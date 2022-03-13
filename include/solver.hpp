#include <math.h> 
#include <vector>
#include <iostream>
#include <execution>
#include <algorithm>
#include <functional>


std::vector<int> indexes;


void start_indexes(int size)
{
    indexes.resize(size);
    for (int i = 0; i < size; i++)
    {
        indexes[i] = i;
    }
    
}


int IX(int x, int y, int N)
{
    return x + y * N;
}



void for_all(int N, std::function<void(int, int)> f)
{
    std::for_each(std::execution::par, indexes.begin(), indexes.end(),
         [&N, &f](int i){
            int y = i / N;
            int x = i % N;
            f(x, y);
    });
}


class Fluid
{
public:
    const int N;
    const int size;
    
    std::vector<double> s;
    std::vector<double> density;

    std::vector<double> vx = {0};
    std::vector<double> vy = {0};

    std::vector<double> vx0 = {0};
    std::vector<double> vy0 = {0};


    Fluid(int N, double dt, double diff, double visc) : 
        N(N), size(N*N), dt(dt), diff(diff), visc(visc)
    {
        start_indexes(this->size);

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

    void step(double dt)
    {
        int N = this->N;
        double visc = this->visc;
        double diff = this->diff;

        this->diffuse(1, this->vx0, this->vx, visc, dt, 4, N);
        this->diffuse(2, this->vy0, this->vy, visc, dt, 4, N);
        
        this->project(this->vx0, this->vy0, this->vx, this->vy, 4, N);
        
        this->advect(1, this->vx, this->vx0, this->vx0, this->vy0, dt, N);
        this->advect(2, this->vy, this->vy0, this->vx0, this->vy0, dt, N);
        
        this->project(this->vx, this->vy, this->vx0, this->vy0, 4, N);
        
        this->diffuse(0, this->s, this->density, diff, dt, 4, N);
        this->advect(0, this->density, this->s, this->vx, this->vy, dt, N);

        this->fade_density();
    }


    void add_density(int x, int y, double density)
    {
        int i = IX(x, y, this->N);
        this->density[i] += density;

        if (std::isnan(this->density[i]))
            this->density[i] = 10;
        if (std::isnan(this->s[i]))
            this->s[i] = 10;

    }


    void add_velocity(int x, int y, double vx, double vy)
    {
        int i = IX(x, y, this->N);

        this->vx[i] += vx;
        this->vy[i] += vy;


        if (std::isnan(this->vx[i]))
        {
            this->vx[i] = 10;
            this->vx0[i] = 10;
        }
            
        if (std::isnan(this->vy[i]))
        {
            this->vy[i] = 10;
            this->vy0[i] = 10;
        }

        if (abs(this->vx[i]) > 10)
            this->vx[i] = this->vx[i] > 0 ? 10 : -10;
        if (abs(this->vy[i]) > 10)
            this->vy[i] = this->vy[i] > 0 ? 10 : -10;

    }


    void clear()
    {
        std::fill(this->vx.begin(), this->vx.end(), 0);
        std::fill(this->vy.begin(), this->vy.end(), 0);
        std::fill(this->vx0.begin(), this->vx0.end(), 0);
        std::fill(this->vy0.begin(), this->vy0.end(), 0);
        std::fill(this->s.begin(), this->s.end(), 0);
        std::fill(this->density.begin(), this->density.end(), 0);
    }

private:
    double dt;
    double diff;
    double visc;


    void lin_solve(int b, std::vector<double>& x, std::vector<double>& x0,
        double a, double c, int iter, int N)
    {
        double c_inv = 1.0 / c;
        
        for (int k=0; k<iter; k++)
        {
            for_all(N, [&a, &N, &x, &x0, &c_inv](int i, int j) {
                double nx = (x0[IX(i, j, N)] + a * (
                        x[IX(i+1, j, N)] + x[IX(i-1, j, N)] +
                        x[IX(i, j+1, N)] + x[IX(i, j-1, N)]
                    )) * c_inv;

                if (abs(nx) < 1000)
                    x[IX(i, j, N)] = nx;
            });
            
            set_boundary(b, x, N);
        }
    }


    void diffuse(int b, std::vector<double>& x, std::vector<double>& x0,
        double diff, double dt, int iter, int N)
    {
        double a = dt * diff * (N - 2) * (N - 2);
        lin_solve(b, x, x0, a, 1 + 4 * a, iter, N);
    }


    void project(std::vector<double>& vx, std::vector<double>& vy, 
        std::vector<double>& p, std::vector<double>& div, int iter, int N)
    {
        for_all(N, [&vx, &vy, &div, &N](int i, int j) {
            div[IX(i, j, N)] = -0.5f*(
                    vx[IX(i+1, j, N)] -vx[IX(i-1, j, N)] +
                    vy[IX(i, j+1, N)] -vy[IX(i, j-1,N)]
            )/N;
        });

        std::fill(p.begin()+1, p.end()-1, 0);

        set_boundary(0, div, N); 
        set_boundary(0, p, N);
        lin_solve(0, p, div, 1, 4, iter, N);
        
        for_all(N, [&](int i, int j){
            vx[IX(i, j, N)] -= 
                0.5f * (p[IX(i+1, j, N)] -p[IX(i-1, j, N)]) * N;
            vy[IX(i, j, N)] -= 
                0.5f * (p[IX(i, j+1, N)] -p[IX(i, j-1, N)]) * N;
        });

        set_boundary(1, vx, N);
        set_boundary(2, vy, N);
    }


    void advect(int b, std::vector<double>& d, std::vector<double>& d0, 
        std::vector<double>& velocX, std::vector<double>& velocY, 
        double dt, int N)
    {
        for_all(N, [&](int i, int j) 
        {
            double x = i - dt * (N - 2) * velocX[IX(i, j, N)];
            double y = j - dt * (N - 2) * velocY[IX(i, j, N)];
            
            auto saturate = [](double v, double l, double h) {
                return std::min(h, std::max(l, v));
            };
            
            int xs = std::floor(saturate(x, 0.5f, N + 0.5f));
            int ys = std::floor(saturate(y, 0.5f, N + 0.5f));
            
            double s1 = x - xs; 
            double s0 = 1.0f - s1; 
            double t1 = y - ys; 
            double t0 = 1.0f - t1;

            d[IX(i, j, N)] = 
                s0 * (t0 * d0[IX(xs, ys, N)] + t1 * d0[IX(xs, ys+1, N)]) +
                s1 * (t0 * d0[IX(xs+1, ys, N)] + t1 * d0[IX(xs+1, ys+1, N)]);
        });

        set_boundary(b, d, N);
    }


    void set_boundary(int b, std::vector<double>& x, int N)
    {
        for(int i = 1; i < N - 1; i++) 
        {
            x[IX(i, 0  , N)] = b == 2 ? -x[IX(i, 1  , N)] : x[IX(i, 1  , N)];
            x[IX(i, N-1, N)] = b == 2 ? -x[IX(i, N-2, N)] : x[IX(i, N-2, N)];
        }
        for(int j = 1; j < N - 1; j++) 
        {
            x[IX(0  , j, N)] = b == 1 ? -x[IX(1  , j, N)] : x[IX(1  , j, N)];
            x[IX(N-1, j, N)] = b == 1 ? -x[IX(N-2, j, N)] : x[IX(N-2, j, N)];
        }
        
        // take the mean value of the boundaries
        x[IX(0, 0, N)] = (x[IX(1, 0, N)] + x[IX(0, 1, N)]) / 2;
        x[IX(0, N-1, N)] = (x[IX(1, N-1, N)] + x[IX(0, N-2, N)]) / 2;
        x[IX(N-1, 0, N)] = (x[IX(N-2, 0, N)] + x[IX(N-1, 1, N)]) / 2;
        x[IX(N-1, N-1, N)] = (x[IX(N-2, N-1, N)] + x[IX(N-1, N-2, N)]) / 2;
    }


    void fade_density()
    {
        for (auto& d:density)
        {
            if (d > 0.1) d-=0.001f;
        }
    }

};