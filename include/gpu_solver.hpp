#include <arrayfire.h>

#include <iostream>

#define IX_0 af::seq(1, af::end-1)


void set_boundary(int b, af::array& x, int N)
{
    x(af::span, 0) = b == 2 ? 0-x(af::span, 1) : x(af::span, 1);
    x(af::span, af::end) = b == 2 ? 0-x(af::span, af::end-1) : x(af::span, af::end-1);

    x(0, af::span) = b == 1 ? 0-x(1, af::span) : x(1, af::span);
    x(af::end, af::span) = b == 1 ? 0-x(af::end-1, af::span) : x(af::end-1, af::span);
    
    x(0, 0) = 0.5f * (x(1, 0) + x(0,1));
    x(0, af::end) = 0.5f * (x(1, af::end) + x(0, af::end-1));
    x(af::end, 0) = 0.5f * (x(af::end-1, 0) + x(af::end, 0));
    x(af::end, af::end) = 0.5f * (x(af::end-1, af::end) + x(af::end, af::end-1));
}


void lin_solve(int b, af::array& x, af::array& x0,
    float a, float c, int iter, int N)
{
    float c_inv = 1.0 / c;
    
    for (int k=0; k<iter; k++)
    {
        x(IX_0, IX_0) = (x0(IX_0, IX_0) + a * (
            x(af::seq(0, af::end-2), IX_0) + x(af::seq(2, af::end), IX_0) +
            x(IX_0, af::seq(0, af::end-2)) + x(IX_0, af::seq(2, af::end))
        ))* c_inv;

        set_boundary(b, x, N);
    }
}


void project(af::array& vx, af::array& vy, 
    af::array& p, af::array& div, int iter, int N)
{
    div(IX_0, IX_0) = -0.5f * (
        vx(af::seq(2, af::end), IX_0) + vx(af::seq(0, af::end-2), IX_0) +
        vy(IX_0, af::seq(2, af::end)) + vy(IX_0, af::seq(0, af::end-2))
    ) / N;
    
    p(IX_0, IX_0) = af::constant(0, p(IX_0, IX_0).dims());

    set_boundary(0, div, N); 
    set_boundary(0, p, N);
    lin_solve(0, p, div, 1, 4, iter, N);

    vx(IX_0, IX_0) -= 0.5f * (
        p(af::seq(2, af::end), IX_0) - p(af::seq(0, af::end-2), IX_0)) * N;
    vy(IX_0, IX_0) -= 0.5f * (
        p(IX_0, af::seq(2, af::end)) - p(IX_0, af::seq(0, af::end-2))) * N;

    set_boundary(1, vx, N);
    set_boundary(2, vy, N);
}


void advect(int b, af::array& d, af::array& d0, 
    af::array& velocX, af::array& velocY, 
    float dt, int N)
{
    float Nfloat = N;
    
    af::array tmp1 = dt * (N - 2) * velocX(IX_0, 0);
    af::array tmp2 = dt * (N - 2) * velocY(IX_0, 0);

    af::array x = af::seq(tmp1.dims(0)) - tmp1;
    af::array y = af::seq(tmp2.dims(0)) - tmp2;

    x = af::min(af::max(0.5, x), Nfloat + 0.5);
    y = af::min(af::max(0.5, y), Nfloat + 0.5);

    af::array i0 = af::floor(x);
    af::array i1 = i0 + 1;

    af::array j0 = floor(y);
    af::array j1 = j0 + 1.0f;

    af::array s1 = af::tile(x - i0, 1, N-2); 
    af::array s0 = 1.0f - s1; 
    af::array t1 = af::tile(y - j0, 1, N-2);
    af::array t0 = 1.0f - t1;
    
    d(IX_0, IX_0) = s0 * (t0 * d0(i0, j0) + t1 * d0(i0, j1)) +
        s1 * (t0 * d0(i1, j0) + t1 * d0(i1, j1));

    set_boundary(b, d, N);
}


class Fluid
{
public:
    const int N;
    
    af::array s;
    af::array density;

    af::array vx;
    af::array vy;

    af::array vx0;
    af::array vy0;


    Fluid(int N, float dt, float diff, float visc) : 
        N(N), dt(dt), diff(diff), visc(visc)
    {
        // init vectors
        this->s = af::constant(0, N, N);
        this->density = af::constant(0, N, N);
        
        this->vx = af::constant(0, N, N);
        this->vy = af::constant(0, N, N);
        this->vx0 = af::constant(0, N, N);
        this->vy0 = af::constant(0, N, N);
    }

    ~Fluid()
    {

    }

    void step()
    {
        int N = this->N;
        float visc = this->visc;
        float diff = this->diff;
        float dt = this->dt;

        diffuse(1, vx0, vx, visc, dt, 4, N);
        diffuse(2, vy0, vy, visc, dt, 4, N);
        
        project(vx0, vy0, vx, vy, 4, N);
        
        advect(1, vx, vx0, vx0, vy0, dt, N);
        advect(2, vy, vy0, vx0, vy0, dt, N);
        
        project(vx, vy, vx0, vy0, 4, N);
        
        diffuse(0, s, density, diff, dt, 4, N);
        advect(0, density, s, vx, vy, dt, N);

        fade_density();
    }


    void add_density(int x, int y, float density)
    {
        this->density(x, y) += density;
    }

    void add_velocity(int x, int y, float vx, float vy)
    {
        this->vx(x, y) += vx;
        this->vy(x, y) += vy;
    }

private:
    float dt;
    float diff;
    float visc;


    void diffuse(int b, af::array& x, af::array& x0,
        float diff, float dt, int iter, int N)
    {
        float a = dt * diff * (N - 2) * (N - 2);
        lin_solve(b, x, x0, a, 1 + 4 * a, iter, N);
    }


    void fade_density()
    {
        this->density = af::max(this->density - 0.0001, 0.001);
    }

};