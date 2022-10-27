#define _WIN32_WINNT 0x0500
#include"ACO.h"
#include <omp.h>
typedef struct
{
    int m1;
    int m2;
    int h1;
    int h2;
    float tx;
    float ty;
    float tz;
    float rx;
    float ry;
    float rz;
    float nCompleted;
    float percentage;
}output;

int main(int argc, char* argv[])
{
    
    //HWND hWnd = GetConsoleWindow();
    //ShowWindow(hWnd, 1);
    flag_cross_section = 3;

    if (argc >= 2)
    {
        //THREADS = (unsigned)argv[1];
        //flag_cross_section = (unsigned)argv[2];
        THREADS = atoi(argv[1]);
        flag_cross_section = atoi(argv[2]);
    }
    else
    {
        THREADS = 8;
        flag_cross_section = 3;
    }
    float tmp_x;
    if (flag_cross_section == choose_arch) tmp_x = Arch::l;
    else if (flag_cross_section == choose_cylinder) tmp_x = Cylinder::r;
    else if (flag_cross_section == choose_horseshoe) tmp_x = Horseshoe::r;
    else tmp_x = Square::r;
    //PSO
    nMid70 = 2;
    nHorizon = 2;
    float p_max[6]{ -(tmp_x - 1.7f), 1.5f, 0.0f,  90.0f,  90.0f,  90.0f };//-5.3
    float p_min[6]{ -(tmp_x - 1.3f), 0.5f, 0.0f, -90.0f, -90.0f, -90.0f };//-5.7
    std::fstream opt("time.csv", std::ios::app | std::ios::out);
    if (!opt.good())std::cerr << "file failed!\n";
    
    for (int i = 11; i >= 11; i--)
    {       
        auto main_time = omp_get_wtime();
        THREADS = i;
        ACO a;
        a.set_p_max(p_max);
        a.set_p_min(p_min);
        a.init();
        a.optimize();
        main_time = omp_get_wtime() - main_time;
        opt << i << "," << main_time << "\n";
        cout << "\nmain() time=" << main_time << "seconds.";
    }
    opt.close();   
    return 0;
}