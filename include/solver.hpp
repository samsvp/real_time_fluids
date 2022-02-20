#include <math.h> 
#include <vector>
#include <iostream>


int IX(int x, int y, int N)
{
    return x + y * N;
}

void set_boundary(int b, std::vector<float>& x, int N)
{
    for(int i = 1; i < N - 1; i++) {
        x[IX(i, 0  , N)] = b == 2 ? -x[IX(i, 1  , N)] : x[IX(i, 1  , N)];
        x[IX(i, N-1, N)] = b == 2 ? -x[IX(i, N-2, N)] : x[IX(i, N-2, N)];
    }
    for(int j = 1; j < N - 1; j++) {
        x[IX(0  , j, N)] = b == 1 ? -x[IX(1  , j, N)] : x[IX(1  , j, N)];
        x[IX(N-1, j, N)] = b == 1 ? -x[IX(N-2, j, N)] : x[IX(N-2, j, N)];
    }
    
    x[IX(0, 0, N)] = 0.5f * (x[IX(1, 0, N)] + x[IX(0, 1, N)]);
    x[IX(0, N-1, N)] = 0.5f * (x[IX(1, N-1, N)] + x[IX(0, N-2, N)]);
    x[IX(N-1, 0, N)] = 0.5f * (x[IX(N-2, 0, N)] + x[IX(N-1, 1, N)]);
    x[IX(N-1, N-1, N)] = 0.5f * (x[IX(N-2, N-1, N)] + x[IX(N-1, N-2, N)]);
}


void lin_solve(int b, std::vector<float>& x, std::vector<float>& x0,
    float a, float c, int iter, int N)
{
    float c_inv = 1.0 / c;
    
    for (int k=0; k<iter; k++)
    {
        for (int j=1; j<N - 1; j++)
        {
            for (int i=1; i<N - 1; i++)
            {
                x[IX(i, j, N)] = (x0[IX(i, j, N)] + a * (
                    x[IX(i+1, j, N)] + x[IX(i-1, j, N)] +
                    x[IX(i, j+1, N)] + x[IX(i, j-1, N)]
                )) * c_inv;
            }
        }
        set_boundary(b, x, N);
    }
}


void project(std::vector<float>& vx, std::vector<float>& vy, 
    std::vector<float>& p, std::vector<float>& div, int iter, int N)
{
    for (int j = 1; j < N - 1; j++) {
        for (int i = 1; i < N - 1; i++) {
            div[IX(i, j, N)] = -0.5f*(
                        vx[IX(i+1, j, N)] -vx[IX(i-1, j, N)] +
                        vy[IX(i, j+1, N)] -vy[IX(i, j-1,N)]
                )/N;
            p[IX(i, j, N)] = 0;
        }
    }

    set_boundary(0, div, N); 
    set_boundary(0, p, N);
    lin_solve(0, p, div, 1, 4, iter, N);
    
    for (int j = 1; j < N - 1; j++) {
        for (int i = 1; i < N - 1; i++) {
            vx[IX(i, j, N)] -= 
                0.5f * (p[IX(i+1, j, N)] -p[IX(i-1, j, N)]) * N;
            vy[IX(i, j, N)] -= 
                0.5f * (p[IX(i, j+1, N)] -p[IX(i, j-1, N)]) * N;
        }
    }

    set_boundary(1, vx, N);
    set_boundary(2, vy, N);
}


void advect(int b, std::vector<float>& d, std::vector<float>& d0, 
    std::vector<float>& velocX, std::vector<float>& velocY, 
    float dt, int N)
{
    float i0, i1, j0, j1;
    
    float dtx = dt * (N - 2);
    float dty = dt * (N - 2);
    
    float s0, s1, t0, t1;
    float tmp1, tmp2, x, y;
    
    float Nfloat = N;
    float ifloat, jfloat;
    int i, j;

    for(j = 1, jfloat = 1; j < N - 1; j++, jfloat++) { 
        for(i = 1, ifloat = 1; i < N - 1; i++, ifloat++) {
            tmp1 = dtx * velocX[IX(i, j, N)];
            tmp2 = dty * velocY[IX(i, j, N)];
            x = ifloat - tmp1; 
            y = jfloat - tmp2;
            
            if(x < 0.5f) x = 0.5f; 
            if(x > Nfloat + 0.5f) x = Nfloat + 0.5f; 
            i0 = floor(x); 
            i1 = i0 + 1.0f;
            if(y < 0.5f) y = 0.5f; 
            if(y > Nfloat + 0.5f) y = Nfloat + 0.5f; 
            j0 = floor(y);
            j1 = j0 + 1.0f;
            
            s1 = x - i0; 
            s0 = 1.0f - s1; 
            t1 = y - j0; 
            t0 = 1.0f - t1;
            
            int i0i = (int)i0;
            int i1i = (int)i1;
            int j0i = (int)j0;
            int j1i = (int)j1;
            

            d[IX(i, j, N)] = 
                s0 * (t0 * d0[IX(i0i, j0i, N)] + t1 * d0[IX(i0i, j1i, N)]) +
                s1 * (t0 * d0[IX(i1i, j0i, N)] + t1 * d0[IX(i1i, j1i, N)]);
        }
    }
    set_boundary(b, d, N);
}


class Fluid
{
public:
    const int N;
    const int size;
    
    std::vector<float> s;
    std::vector<float> density;

    std::vector<float> vx = {0};
    std::vector<float> vy = {0};

    std::vector<float> vx0 = {0};
    std::vector<float> vy0 = {0};


    Fluid(int N, float dt, float diff, float visc) : 
        N(N), size(N*N), dt(dt), diff(diff), visc(visc)
    {
        // init vectors
        this->s.resize(this->size);
        this->density.resize(this->size);
        
        this->vx.resize(this->size);
        this->vy.resize(this->size);
        this->vx0.resize(this->size);
        this->vy0.resize(this->size);
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
        int i = IX(x, y, this->N);
        this->density[i] += density;
    }

    void add_velocity(int x, int y, float vx, float vy)
    {
        int i = IX(x, y, this->N);

        this->vx[i] += vx;
        this->vy[i] += vy;
    }

private:
    float dt;
    float diff;
    float visc;


    void diffuse(int b, std::vector<float>& x, std::vector<float>& x0,
        float diff, float dt, int iter, int N)
    {
        float a = dt * diff * (N - 2) * (N - 2);
        lin_solve(b, x, x0, a, 1 + 4 * a, iter, N);
    }


    void fade_density()
    {
        for (auto& d:density)
        {
            if (d > 1) d-=0.001f;
        }
    }

};