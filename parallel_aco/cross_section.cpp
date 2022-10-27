#ifdef VTK
#include <vtkAutoInit.h>
VTK_MODULE_INIT(vtkRenderingOpenGL2);
VTK_MODULE_INIT(vtkInteractionStyle);
#include <vtkActor.h>
#include <vtkCamera.h>
#include <vtkCellArray.h>
#include <vtkFloatArray.h>
#include <vtkPointData.h>
#include <vtkPoints.h>
#include <vtkPolyData.h>
#include <vtkPolyDataMapper.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkRenderer.h>
#include <vtkInteractorStyleTrackballActor.h>
#include <vtkExtractEdges.h>
#endif
#include "cross_section.h"
#include<iostream>
//#include <math.h>
using namespace std;

float t, t_s, t_p;
pcl::PointCloud<pcl::PointXYZ> cloud;
pcl::PointCloud<pcl::PointXYZ> cloud_cube;

Cross_Section::Cross_Section() {}
std::vector<unsigned>::iterator Cross_Section::check_continue(std::vector<unsigned>::iterator i, std::vector<unsigned>& vec)
{
    if (i + 1 != vec.end() && *(i + 1) - *i != 1) return i;
    else if (i + 1 != vec.end()) return check_continue(i + 1, vec);
    else return i;
}
unsigned Cross_Section::cube_x = 0U;
unsigned Cross_Section::cube_y = 0U;
unsigned Cross_Section::cube_z = 0U;
unsigned Cross_Section::nCube = 0U;
unsigned* Cross_Section::cube = nullptr;

Square::Square() {}

unsigned Square::square(std::vector<Eigen::Vector3d>& v, unsigned length) {
    cloud.width = length - 1 + cloud.size();
    cloud.height = 1;
    cloud.reserve(cloud.width * cloud.height);

    float ox = v[0].x(), oy = v[0].y(), oz = v[0].z();
    unsigned intersection = 0;
    for (unsigned i = 1; i < length; i++)
    {
        float tmpX = 0.0f, tmpY = 0.0f, tmpZ = 0.0f;
        float bx = v[i].x(), by = v[i].y(), bz = v[i].z();
        //求斜率
        float D = by - oy;        
        if (D != 0) {
            float t1 = 0 - oy / D;//y=0
            float tmpX1 = ox + (bx - ox) * t1;
            float tmpY1 = oy + (by - oy) * t1;
            float tmpZ1 = oz + (bz - oz) * t1;
            if (tmpX1 >= -r && tmpX1 <= r)
            {
                tmpX = tmpX1;
                tmpY = 0;
                tmpZ = tmpZ1;              
                pcl::PointXYZ tmp{ tmpX,tmpY,tmpZ };
                cloud.emplace_back(tmp);
                intersection++;
            }            
            float t2 = (l - oy) / D;//y=l
            float tmpX2 = ox + (bx - ox) * t2;
            float tmpY2 = oy + (by - oy) * t2;
            float tmpZ2 = oz + (bz - oz) * t2;
            if (tmpX2 >= -r && tmpX2 <= r)
            {
                tmpX = tmpX2;
                tmpY = l;
                tmpZ = tmpZ2;
                pcl::PointXYZ tmp{ tmpX,tmpY,tmpZ };
                cloud.emplace_back(tmp);
                intersection++;
            }
            //cout << "\ntmpx,y,z=" << tmpX << ", " << tmpY << ", " << tmpZ;
        }
        D = bx - ox;
        if (D != 0) {
            float t1 = (r - ox) / D;//x=r
            float tmpX1 = ox + (bx - ox) * t1;
            float tmpY1 = oy + (by - oy) * t1;
            float tmpZ1 = oz + (bz - oz) * t1;
            if (tmpY1 >= 0 && tmpY1 <= l) {
                tmpX = r;
                tmpY = tmpY1;
                tmpZ = tmpZ1;
                pcl::PointXYZ tmp{ tmpX,tmpY,tmpZ };
                cloud.emplace_back(tmp);
                intersection++;
            }            
            float t2 = (-r - ox) / D;//x=-r
            float tmpX2 = ox + (bx - ox) * t2;
            float tmpY2 = oy + (by - oy) * t2;
            float tmpZ2 = oz + (bz - oz) * t2;
            if (tmpY2 >= 0 && tmpY2 <= l) {
                tmpX = -r;
                tmpY = tmpY2;
                tmpZ = tmpZ2;
                pcl::PointXYZ tmp{ tmpX,tmpY,tmpZ };
                cloud.emplace_back(tmp);
                intersection++;
            }
            //cout << "\ntmpx,y,z=" << tmpX << ", " << tmpY << ", " << tmpZ;
        }      
    }
    // cout << "\n" << __func__ << ":" << intersection << ", cloud-size: " << cloud.size();
    return intersection;
}

float Square::statistics(bool print)
{
    cloud_cube.clear();
    Cross_Section::cube_x = (unsigned)(r * 2.0f * scale) + 1;//51
    Cross_Section::cube_y = (unsigned)(l * scale) + 1;
    Cross_Section::cube_z = (unsigned)(seg * scale) + 1;
    Cross_Section::nCube = (cube_x) * (cube_y) * (cube_z);//50*50*50
    cloud_cube.reserve(nCube);
    Cross_Section::cube = new unsigned[nCube];
    unsigned cloud_size = cloud.size();
    if (print)
    {
        std::cout << "------------------------------------------------------------------------------\n"
            << "Number of points = " << cloud_size << endl;
        std::cout << "Number of cubes = " << cube_x << " * " << cube_y << " * " << cube_z << " = " << nCube
            << "\nSegmentation: " << std::fixed << setprecision(2) << end - start << "m, precision: " << precision << "m\n";
    }
    #pragma omp parallel for num_threads(THREADS) schedule(dynamic, nCube/THREADS)
    for (unsigned i = 0; i < nCube; i++)
    {
        cube[i] = 0U;
    }
    auto start_1 = std::chrono::steady_clock::now();
    #pragma omp parallel for num_threads(THREADS) schedule(dynamic, cloud_size/THREADS)
    for (unsigned i = 0; i < cloud_size; i++)
    {
        if (cloud.points.at(i).z <= end && cloud.points.at(i).z >= start)
        {
            unsigned _x = (unsigned)((cloud.points.at(i).x + r) * scale);//-2.5->0  //2.5->50
            unsigned _y = (unsigned)((cloud.points.at(i).y) * scale);//0
            unsigned _z = (unsigned)((cloud.points.at(i).z - start) * scale);//0
            ++cube[_z * cube_x * cube_y+_x * cube_y+_y];//0 //1*51 //2*51
            //cout << "\nx,y,z= " << cloud.points.at(i).x << "," << cloud.points.at(i).y << "," << cloud.points.at(i).z;
            //cout << "\n_x,_y,_z= " << _x << "," << _y << "," << _z;
        }
    }
    auto end_1 = std::chrono::duration<float>(std::chrono::steady_clock::now() - start_1);
    t_p += end_1.count();

    unsigned long long sCube = 0, sPts = 0, sN = 0, sVar = 0;

    float threshold = precision * sqrt(2.0) * 0.9f;
    unsigned startScale = start * scale, endScale = end * scale;
    unsigned nCompleted = 0, countCompleted = 0;
    std::vector<unsigned> segCompleted, vec_max_nPts, vec_min_nPts;
    std::vector<float> vec_ring_ave, vec_ring_var;
    vec_max_nPts.reserve(cube_z);
    vec_min_nPts.reserve(cube_z);
    vec_ring_ave.reserve(cube_z);
    vec_ring_var.reserve(cube_z);
    segCompleted.reserve(cube_z);

    for (unsigned k = 0; k < cube_z; k++)
    {
        float z_ = (float)(k + start * scale) * precision;
        unsigned completed = 0;
        unsigned cross_section_cube = 0;
        unsigned max_nPts = 0;
        unsigned min_nPts = UINT_MAX;
        float ringAverage = 0.0f, ringVar = 0.0f;
        for (unsigned i = 0; i < cube_x; i++)//0-50 個方塊
        {
            float x_ = (float)(i - r * scale) * precision;//(50-2.5*10)*0.1 = 2.5 //(0-2.5*10)*0.1=-2.5
            for (unsigned j = 0; j < cube_y; j++)
            {
                float y_ = (float)(j)*precision;//50 * 0.1 = 5
                //cout << "\ny_=" << y_;
                if (x_ <= r && x_ >= -r)//-2.5 < x < 2.5
                {
                    if (abs(y_-l) <= eps||y_ < eps)
                    {
                        if (print) {
                            pcl::PointXYZ tmp{ x_,y_,z_ };
                            cloud_cube.emplace_back(tmp);
                        }
                        ++sCube;
                        ++cross_section_cube;
                        auto n = (unsigned long long)cube[k * (cube_x) * (cube_y)+i * (cube_y)+j];
                        //opt << n << ", ";
                        if (n > max_nPts)max_nPts = n;
                        if (n < min_nPts)min_nPts = n;
                        if (n)
                        {
                            if (z_ >= 5.0f && z_ <= 15.0f)
                            {
#ifdef VTK
                                vtkPolyData* v_cube = vtkPolyData::New();
                                vtkCellArray* polys = vtkCellArray::New();
                                vtkPoints* points = vtkPoints::New();
                                vtkPolyDataMapper* cubeMapper = vtkPolyDataMapper::New();
                                cubeMapper->SetScalarRange(0, 1);
                                for (int p = 0; p < 8; p++)
                                {
                                    points->InsertPoint(p, corner[p]);
                                }
                                for (int p = 0; p < 6; p++)
                                {
                                    polys->InsertNextCell(4, pts[p]);
                                }
                                v_cube->SetPoints(points);
                                v_cube->SetPolys(polys);
                                vtkFloatArray* scalars = vtkFloatArray::New();
                                for (int p = 0; p < 8; p++)
                                {
                                    scalars->InsertTuple1(p, (float)n / 10000.0f);
                                }
                                v_cube->GetPointData()->SetScalars(scalars);
                                //scalars->Delete();

                                //vtkExtractEdges* extract = vtkExtractEdges::New();
                                //extract->SetInputData(cube);
                                //vtkPolyDataMapper* mapEdges = vtkPolyDataMapper::New();
                                //mapEdges->SetInputConnection(extract->GetOutputPort());
                                //mapEdges->SetScalarVisibility(0);
                                //vtkActor* edgeActor = vtkActor::New();
                                //edgeActor->SetMapper(mapEdges);
                                //edgeActor->VisibilityOn();

                                cubeMapper->SetInputData(v_cube);

                                vtkActor* cubeActor = vtkActor::New();
                                cubeActor->SetPosition(x_, y_, z_);
                                vtkActor* tempactor = vtkActor::New();
                                cubeActor->SetMapper(cubeMapper);
                                renderer->AddActor(cubeActor);
#endif
                            }
                            ringVar += n * n;
                            sVar += n * n;
                            ringAverage += n;
                            sPts += n;
                            ++sN;
                            ++completed;
                        }
                    }
                }
                if (y_ >= 0 && y_ <= l)
                {
                    if (abs(x_ - r) <= eps||abs(x_ + r) <= eps)//x=r
                    {
                        if (print) {
                            pcl::PointXYZ tmp{ x_,y_,z_ };
                            cloud_cube.emplace_back(tmp);
                        }
                        ++sCube;
                        ++cross_section_cube;
                        auto n = (unsigned long long)cube[k * (cube_x) * (cube_y)+i * (cube_y)+j];
                        //opt << n << ", ";
                        if (n > max_nPts)max_nPts = n;
                        if (n < min_nPts)min_nPts = n;
                        if (n)
                        {
                            if (z_ >= 5.0f && z_ <= 15.0f)
                            {
#ifdef VTK
                                vtkPolyData* v_cube = vtkPolyData::New();
                                vtkCellArray* polys = vtkCellArray::New();
                                vtkPoints* points = vtkPoints::New();
                                vtkPolyDataMapper* cubeMapper = vtkPolyDataMapper::New();
                                cubeMapper->SetScalarRange(0, 1);
                                for (int p = 0; p < 8; p++)
                                {
                                    points->InsertPoint(p, corner[p]);
                                }
                                for (int p = 0; p < 6; p++)
                                {
                                    polys->InsertNextCell(4, pts[p]);
                                }
                                v_cube->SetPoints(points);
                                v_cube->SetPolys(polys);
                                vtkFloatArray* scalars = vtkFloatArray::New();
                                for (int p = 0; p < 8; p++)
                                {
                                    scalars->InsertTuple1(p, (float)n / 10000.0f);
                                }
                                v_cube->GetPointData()->SetScalars(scalars);
                                //scalars->Delete();

                                //vtkExtractEdges* extract = vtkExtractEdges::New();
                                //extract->SetInputData(cube);
                                //vtkPolyDataMapper* mapEdges = vtkPolyDataMapper::New();
                                //mapEdges->SetInputConnection(extract->GetOutputPort());
                                //mapEdges->SetScalarVisibility(0);
                                //vtkActor* edgeActor = vtkActor::New();
                                //edgeActor->SetMapper(mapEdges);
                                //edgeActor->VisibilityOn();

                                cubeMapper->SetInputData(v_cube);

                                vtkActor* cubeActor = vtkActor::New();
                                cubeActor->SetPosition(x_, y_, z_);
                                vtkActor* tempactor = vtkActor::New();
                                cubeActor->SetMapper(cubeMapper);
                                renderer->AddActor(cubeActor);
#endif
                            }
                            ringVar += n * n;
                            sVar += n * n;
                            ringAverage += n;
                            sPts += n;
                            ++sN;
                            ++completed;
                        }
                    }
                }
                //                if (y_ == 5)
                //                {
                //                    if (print) {
                //                        pcl::PointXYZ tmp{ x_,y_,z_ };
                //                        cloud_cube.emplace_back(tmp);
                //                    }
                //                    ++sCube;
                //                    ++cross_section_cube;
                //                    auto n = (unsigned long long)cube[k * (cube_x) * (cube_y) + i * (cube_y) + j];
                //                    //opt << n << ", ";
                //                    if (n > max_nPts)max_nPts = n;
                //                    if (n < min_nPts)min_nPts = n;
                //                    if (n)
                //                    {
                //                        if (z_ >= 5.0f && z_ <= 15.0f)
                //                        {
                //#ifdef VTK
                //                            vtkPolyData* v_cube = vtkPolyData::New();
                //                            vtkCellArray* polys = vtkCellArray::New();
                //                            vtkPoints* points = vtkPoints::New();
                //                            vtkPolyDataMapper* cubeMapper = vtkPolyDataMapper::New();
                //                            cubeMapper->SetScalarRange(0, 1);
                //                            for (int p = 0; p < 8; p++)
                //                            {
                //                                points->InsertPoint(p, corner[p]);
                //                            }
                //                            for (int p = 0; p < 6; p++)
                //                            {
                //                                polys->InsertNextCell(4, pts[p]);
                //                            }
                //                            v_cube->SetPoints(points);
                //                            v_cube->SetPolys(polys);
                //                            vtkFloatArray* scalars = vtkFloatArray::New();
                //                            for (int p = 0; p < 8; p++)
                //                            {
                //                                scalars->InsertTuple1(p, (float)n / 10000.0f);
                //                            }
                //                            v_cube->GetPointData()->SetScalars(scalars);
                //                            //scalars->Delete();
                //
                //                            //vtkExtractEdges* extract = vtkExtractEdges::New();
                //                            //extract->SetInputData(cube);
                //                            //vtkPolyDataMapper* mapEdges = vtkPolyDataMapper::New();
                //                            //mapEdges->SetInputConnection(extract->GetOutputPort());
                //                            //mapEdges->SetScalarVisibility(0);
                //                            //vtkActor* edgeActor = vtkActor::New();
                //                            //edgeActor->SetMapper(mapEdges);
                //                            //edgeActor->VisibilityOn();
                //
                //                            cubeMapper->SetInputData(v_cube);
                //
                //                            vtkActor* cubeActor = vtkActor::New();
                //                            cubeActor->SetPosition(x_, y_, z_);
                //                            vtkActor* tempactor = vtkActor::New();
                //                            cubeActor->SetMapper(cubeMapper);
                //                            renderer->AddActor(cubeActor);
                //#endif
                //                        }
                //                        ringVar += n * n;
                //                        sVar += n * n;
                //                        ringAverage += n;
                //                        sPts += n;
                //                        ++sN;
                //                        ++completed;
                //                    }
                //                }
                //                else if (y_ <= eps)
                //                {
                //                    if (print) {
                //                        pcl::PointXYZ tmp{ x_,y_,z_ };
                //                        cloud_cube.emplace_back(tmp);
                //                    }
                //                    ++sCube;
                //                    ++cross_section_cube;
                //                    auto n = (unsigned long long)cube[k * cube_x * cube_y + i * cube_y + j];
                //                    //opt << n << ", ";
                //                    if (n)
                //                    {
                //                        if (z_ >= 5.0f && z_ <= 15.0f)
                //                        {
                //#ifdef VTK
                //                            vtkPolyData* v_cube = vtkPolyData::New();
                //                            vtkCellArray* polys = vtkCellArray::New();
                //                            vtkPoints* points = vtkPoints::New();
                //                            vtkPolyDataMapper* cubeMapper = vtkPolyDataMapper::New();
                //                            cubeMapper->SetScalarRange(0, 1);
                //                            for (int p = 0; p < 8; p++)
                //                            {
                //                                points->InsertPoint(p, corner[p]);
                //                            }
                //                            for (int p = 0; p < 6; p++)
                //                            {
                //                                polys->InsertNextCell(4, pts[p]);
                //                            }
                //                            v_cube->SetPoints(points);
                //                            v_cube->SetPolys(polys);
                //                            vtkFloatArray* scalars = vtkFloatArray::New();
                //                            for (int p = 0; p < 8; p++)
                //                            {
                //                                scalars->InsertTuple1(p, (float)n / 10000.0f);
                //                            }
                //                            v_cube->GetPointData()->SetScalars(scalars);
                //                            //scalars->Delete();
                //
                //                            //vtkExtractEdges* extract = vtkExtractEdges::New();
                //                            //extract->SetInputData(cube);
                //                            //vtkPolyDataMapper* mapEdges = vtkPolyDataMapper::New();
                //                            //mapEdges->SetInputConnection(extract->GetOutputPort());
                //                            //mapEdges->SetScalarVisibility(0);
                //                            //vtkActor* edgeActor = vtkActor::New();
                //                            //edgeActor->SetMapper(mapEdges);
                //                            //edgeActor->VisibilityOn();
                //
                //                            cubeMapper->SetInputData(v_cube);
                //
                //                            vtkActor* cubeActor = vtkActor::New();
                //                            cubeActor->SetPosition(x_, y_, z_);
                //                            vtkActor* tempactor = vtkActor::New();
                //                            cubeActor->SetMapper(cubeMapper);
                //                            renderer->AddActor(cubeActor);
                //#endif
                //                        }
                //                        sVar += n * n;
                //                        sPts += n;
                //                        ++sN;
                //                        ++completed;
                //                    }
                //                    //continue;
                //                }
                //                else if (abs(x_) >= r - eps && y_ <=5 )
                //                {
                //                    if (print) {
                //                        pcl::PointXYZ tmp{ x_,y_,z_ };
                //                        cloud_cube.emplace_back(tmp);
                //                    }
                //                    ++sCube;
                //                    ++cross_section_cube;
                //                    auto n = (unsigned long long)cube[k * cube_x * cube_y + i * cube_y + j];
                //                    if (n)
                //                    {
                //                        if (z_ >= 5.0f && z_ <= 15.0f)
                //                        {
                //#ifdef VTK
                //                            vtkPolyData* v_cube = vtkPolyData::New();
                //                            vtkCellArray* polys = vtkCellArray::New();
                //                            vtkPoints* points = vtkPoints::New();
                //                            vtkPolyDataMapper* cubeMapper = vtkPolyDataMapper::New();
                //                            cubeMapper->SetScalarRange(0, 1);
                //                            for (int p = 0; p < 8; p++)
                //                            {
                //                                points->InsertPoint(p, corner[p]);
                //                            }
                //                            for (int p = 0; p < 6; p++)
                //                            {
                //                                polys->InsertNextCell(4, pts[p]);
                //                            }
                //                            v_cube->SetPoints(points);
                //                            v_cube->SetPolys(polys);
                //                            vtkFloatArray* scalars = vtkFloatArray::New();
                //                            for (int p = 0; p < 8; p++)
                //                            {
                //                                scalars->InsertTuple1(p, (float)n / 10000.0f);
                //                            }
                //                            v_cube->GetPointData()->SetScalars(scalars);
                //                            //scalars->Delete();
                //
                //                            //vtkExtractEdges* extract = vtkExtractEdges::New();
                //                            //extract->SetInputData(cube);
                //                            //vtkPolyDataMapper* mapEdges = vtkPolyDataMapper::New();
                //                            //mapEdges->SetInputConnection(extract->GetOutputPort());
                //                            //mapEdges->SetScalarVisibility(0);
                //                            //vtkActor* edgeActor = vtkActor::New();
                //                            //edgeActor->SetMapper(mapEdges);
                //                            //edgeActor->VisibilityOn();
                //
                //                            cubeMapper->SetInputData(v_cube);
                //
                //                            vtkActor* cubeActor = vtkActor::New();
                //                            cubeActor->SetPosition(x_, y_, z_);
                //                            vtkActor* tempactor = vtkActor::New();
                //                            cubeActor->SetMapper(cubeMapper);
                //                            renderer->AddActor(cubeActor);
                //#endif
                //                        }
                //                        sVar += n * n;
                //                        sPts += n;
                //                        ++sN;
                //                        ++completed;
                //                    }
                //                    //continue;
                //                }
            }
        }
        //opt << "\n";
        ringAverage /= cross_section_cube;
        float ringVariance = sqrt(ringVar / (float)cross_section_cube - ringAverage * ringAverage);
        vec_ring_ave.emplace_back(ringAverage);
        vec_ring_var.emplace_back(ringVariance);
        vec_max_nPts.emplace_back(max_nPts);
        vec_min_nPts.emplace_back(min_nPts);
        cout <<"\n" << completed << "/" << cross_section_cube ;
        if ((float)completed / (float)cross_section_cube >= 0.9f)
        {
            ++countCompleted;
            segCompleted.emplace_back((unsigned)(z_ * scale));
        }
        else
        {
            countCompleted = 0U;
        }
        if (countCompleted > nCompleted)nCompleted = countCompleted;
    }

#ifdef VTK
    vtkCamera* camera = vtkCamera::New();
    camera->SetPosition(0, 0, -10);
    camera->SetFocalPoint(0, 0, 0);


    vtkRenderWindow* reWin = vtkRenderWindow::New();
    reWin->AddRenderer(renderer);

    vtkRenderWindowInteractor* iren = vtkRenderWindowInteractor::New();
    iren->SetRenderWindow(reWin);
    vtkInteractorStyleTrackballActor* style = vtkInteractorStyleTrackballActor::New();
    iren->SetInteractorStyle(style);

    renderer->SetActiveCamera(camera);
    renderer->ResetCamera();
    renderer->SetBackground(1, 1, 1);

    reWin->SetSize(300, 300);

    reWin->Render();
    iren->Initialize();
    iren->Start();

    points->Delete();
    v_cube->Delete();
    cubeMapper->Delete();
    renderer->Delete();
    reWin->Delete();
    iren->Delete();
    polys->Delete();
#endif
    //opt.close();
    if (segCompleted.size() == 0)segCompleted.emplace_back(0U);

    if (print)
    {
        //這不應該在這裡吧，會害後面分不同 Sensor 顏色輸出時壞掉...
        //pcl::PassThrough<pcl::PointXYZ> pass;
        //pass.setInputCloud(cloud.makeShared());
        //pass.setFilterFieldName("z");
        //pass.setFilterLimits(0.0, 50.0);
        //pass.filter(cloud);
        //pcl::io::savePCDFileBinaryCompressed<pcl::PointXYZ>("simulated_cylinder.pcd", cloud);
        float sAverage = (float)sPts / (float)sCube;
        float sStandard = sqrt((float)sVar / (float)sCube - sAverage * sAverage);
        std::cout << "------------------------------------------------------------------------------\nStatistics\n";
        std::cout << "\nSum: " << sPts << "(" << std::fixed << setprecision(2) << (float)sPts / cloud_size * 100.0f
            << "%), " << sPts / (r * (2 + EIGEN_PI) * seg * 10000.0f) << "pts/cm^2\n"
            << "Cubes: " << sCube << "(Zero: " << sCube - sN << "(" << (float)(sCube - sN) / (float)sCube * 100.0f << "%), N: "
            << sN << "(" << (float)sN / (float)sCube * 100.0f << "%))";
        std::cout << "\n g = " << std::fixed << setprecision(3) << sAverage << ",  m = " << sStandard;
        std::cout << "\nnCompleted = " << std::fixed << setprecision(1) << nCompleted * precision << endl;
        for (auto i = segCompleted.begin(); i != segCompleted.end();)
        {
            auto it = check_continue(i, segCompleted);
            std::cout << setw(4) << *i * precision << " - " << setw(4) << *it * precision << "\n";
            i = it + 1;
        }
        //std::cout << "------------------------------------------------------------------------------\n";
        //std::cout << setw(6) << "z" << setw(10) << "ave" << setw(10) << "var" << setw(10) << "max" << setw(10) << "min" << endl;
        //for (size_t i = 0; i < vec_max_nPts.size(); i++)
        //{
        //  std::cout << setw(5) << (float)i * precision << "m" << setw(10) << vec_ring_ave.at(i) << setw(10) << vec_ring_var.at(i) << setw(10) << vec_max_nPts.at(i) << setw(10) << vec_min_nPts.at(i) << "\n";
        //}
        //std::cout << "------------------------------------------------------------------------------\n\n";

        //std::fstream opt("diff_amount.csv", std::ios::app | std::ios::out);
        //if (!opt.good())std::cerr << "file failed!\n";
        //opt << "\n, Sum, " << sPts << ", " << (float)sPts / cloud_size * 100.0f << "%, "
        //  << sPts / (r * (2.0f + EIGEN_PI) * seg * 10000.0f) << "pts/cm^2, \n"
        //  << ", Cubes, " << sCube << ", Zero, " << sCube - sN << ", " << (float)(sCube - sN) / sCube * 100.0f << "%, N, "
        //  << sN << ", " << (float)sN / sCube * 100.0f << "%, \n"
        //  << ",  g, " << sAverage << ",  m, " << sStandard << ", \n"
        //  << ", nCompleted, " << nCompleted * precision << *segCompleted.cbegin() * precision << ", \n";
        //opt.close();
    }

    delete[] cube;
    cout << "\n# points: " << cloud.size();
    return (float)nCompleted * precision;
}

Cylinder::Cylinder() {}

unsigned Cylinder::cylinder(std::vector<Eigen::Vector3d>& v, unsigned length) {
    cloud.width = length - 1 + cloud.size();
    cloud.height = 1;
    cloud.reserve(cloud.width * cloud.height);

    float ox = v[0].x(), oy = v[0].y(), oz = v[0].z();
    unsigned intersection = 0;
    for (unsigned i = 0; i < length; i++)
    {
        float tmpX = 0.0f, tmpY = 0.0f, tmpZ = 0.0f;
        float bx = v[i].x(), by = v[i].y(), bz = v[i].z();
        //Ax+By+C, x^2+y^2=r^2
        float A = (bx - ox) * (bx - ox) + (by - oy) * (by - oy);
        float B = 2.0f * (ox * (bx - ox) + oy * (by - oy));
        float C = ox * ox + oy * oy - r * r;
        float D = B * B - 4 * A * C;
        if (D > 0.0f)
        {
            float t1 = (-B + sqrt(D)) / (A * 2.0f)/*, t2 = (-B - sqrt(D)) / (A * 2.0f)*/;
            float tmpX1 = ox + (bx - ox) * t1/*, tmpX2 = ox + (bx - ox) * t2*/;
            float tmpY1 = oy + (by - oy) * t1/*, tmpY2 = oy + (by - oy) * t2*/;
            float tmpZ1 = oz + (bz - oz) * t1/*, tmpZ2 = oz + (bz - oz) * t2*/;
            if (tmpY1 >= h)
            {
                tmpX = tmpX1;
                tmpY = tmpY1;
                tmpZ = tmpZ1;
            }
            else if (tmpY1 - h < eps)
            {
                //float t = oy / (oy - by);
                float t = (h - oy) / (by - oy);
                tmpX = ox + (bx - ox) * t;
                tmpY = h;
                tmpZ = oz + (bz - oz) * t;
            }
            pcl::PointXYZ tmp{ tmpX,tmpY,tmpZ };
            cloud.emplace_back(tmp);
            intersection++;
        }
    }
    // cout << "\n" << __func__ << ":" << intersection << ", cloud-size: " << cloud.size();
    return intersection;
}

float Cylinder::statistics(bool print)//統計
{
    cloud_cube.clear();//放點的空間
    Cross_Section::cube_x = (unsigned)(r * 2.0f * scale) + 1U;//x方向的方塊 //種樹問題
    Cross_Section::cube_y = (unsigned)((r + h) * scale) + 1U;
    Cross_Section::cube_z = (unsigned)(seg * scale) + 1U;
    Cross_Section::nCube = cube_x * cube_y * cube_z;
    cloud_cube.reserve(nCube);
    Cross_Section::cube = new unsigned[nCube];//放方塊的空間

    unsigned cloud_size = cloud.size();//有幾個點
    if (print)
    {
        std::cout << "------------------------------------------------------------------------------\n"
            << "Number of points = " << cloud_size << endl;
        std::cout << "Number of cubes = " << cube_x << " * " << cube_y << " * " << cube_z << " = " << nCube
            << "\nSegmentation: " << std::fixed << setprecision(2) << end - start << "m, precision: " << precision << "m\n";
    }

    #pragma omp parallel for num_threads(THREADS) schedule(dynamic, nCube/THREADS)
    for (unsigned i = 0; i < nCube; i++)
    {
        cube[i] = 0U;
    }//每個放方塊的空間都填0

    auto start_1 = std::chrono::steady_clock::now();
    #pragma omp parallel for num_threads(THREADS) schedule(dynamic, cloud_size/THREADS)
    for (unsigned i = 0; i < cloud_size; i++)
    {
        if (cloud.points.at(i).z <= end && cloud.points.at(i).z >= start)//如果在位置i上的z在範圍內
        {
            unsigned _x = (unsigned)((cloud.points.at(i).x + r) * scale);//方塊的x索引值
            unsigned _y = (unsigned)((cloud.points.at(i).y + h) * scale);//從原點到y有多少個方塊
            unsigned _z = (unsigned)((cloud.points.at(i).z - start) * scale);//從原點到z有多少個方塊
            //cout << "\n_x=" << _x << "   _y=" << _y << "   _z=" << _z;
            ++cube[_z * cube_x * cube_y + _x * cube_y + _y];//為了數多少個點打在同個方塊上
        }
    }
    auto end_1 = std::chrono::duration<float>(std::chrono::steady_clock::now() - start_1);
    t_p += end_1.count();


    unsigned long long sCube = 0, sPts = 0, sN = 0, sVar = 0;

    float threshold = precision * sqrt(2.0) * 0.9f;//臨界點
    unsigned startScale = start * scale, endScale = end * scale;//50公尺*一公尺有幾個方塊
    unsigned nCompleted = 0, countCompleted = 0;
    std::vector<unsigned> segCompleted, vec_max_nPts, vec_min_nPts;
    std::vector<float> vec_ring_ave, vec_ring_var;
    vec_max_nPts.reserve(cube_z);
    vec_min_nPts.reserve(cube_z);
    vec_ring_ave.reserve(cube_z);
    vec_ring_var.reserve(cube_z);
    segCompleted.reserve(cube_z);
    //std::fstream opt("var.csv", std::ios::app | std::ios::out);
    //vtkRenderer* renderer = vtkRenderer::New();
    //float corner[8][3] = { { 0, 0, 0 }, { precision, 0, 0 }, { precision, 0, precision }, { 0, 0, precision }, { 0, precision, 0 }, { precision, precision, 0 }, { precision, precision, precision }, { 0, precision, precision } };
    //vtkIdType pts[6][4] = { { 0, 1, 2, 3 }, { 4, 5, 6, 7 }, { 2, 3, 7, 6 }, { 1, 0, 4, 5 }, { 0, 4, 7, 3 }, { 1, 2, 6, 5 } };

    for (unsigned k = 0; k < cube_z; k++)
    {
        float z_ = (float)(k + start * scale) * precision;//0.1公尺一方塊 //方塊的距離
        unsigned completed = 0;
        unsigned cross_section_cube = 0;//形狀分成的方塊
        unsigned max_nPts = 0;
        unsigned min_nPts = UINT_MAX;
        float ringAverage = 0.0f, ringVar = 0.0f;
        for (unsigned i = 0; i < cube_x; i++)
        {
            float x_ = (float)(i - r * scale) * precision;//方塊的x距離
            for (unsigned j = 0; j < cube_y; j++)
            {
                float y_ = (float)(j - h * scale) * precision;
                float r_after = (x_ < 0.0f)
                    ? sqrt((x_ + precision) * (x_ + precision) + y_ * y_)
                    : sqrt(x_ * x_ + y_ * y_);
                if ((int&)y_ <= (int&)h || (((int&)r >= (int&)r_after) && abs(r - r_after) - threshold <= eps))
                {//y<0 或 z<r 或 閥值=不合理 就不做
                    if (print) {
                        pcl::PointXYZ tmp{ x_,y_,z_ };
                        cloud_cube.emplace_back(tmp);
                    }
                    ++sCube;//被點雲覆蓋的方塊++
                    ++cross_section_cube;//組成完整環的方塊++
                    auto n = (unsigned long long)cube[k * cube_x * cube_y + i * cube_y + j];
                    //opt << n << ", ";
                    if (n > max_nPts)max_nPts = n;
                    if (n < min_nPts)min_nPts = n;
                    if (n)
                    {
                        ringVar += n * n;
                        sVar += n * n;
                        ringAverage += n;
                        sPts += n;
                        ++sN;
                        ++completed;
                    }
                }
            }
        }
        //opt << "\n";
        ringAverage /= cross_section_cube;
        float ringVariance = sqrt(ringVar / (float)cross_section_cube - ringAverage * ringAverage);
        vec_ring_ave.emplace_back(ringAverage);
        vec_ring_var.emplace_back(ringVariance);
        vec_max_nPts.emplace_back(max_nPts);
        vec_min_nPts.emplace_back(min_nPts);
        cout << completed << "/" << cross_section_cube << endl;
        if ((float)completed / (float)cross_section_cube >= 0.9f)
        {
            ++countCompleted;//完整環++
            segCompleted.emplace_back((unsigned)(z_ * scale));//最長連續完整環長度
        }
        else
        {
            countCompleted = 0U;
        }
        if (countCompleted > nCompleted)nCompleted = countCompleted;
    }
#ifdef VTK
    vtkCamera* camera = vtkCamera::New();
    camera->SetPosition(0, 0, -10);
    camera->SetFocalPoint(0, 0, 0);


    vtkRenderWindow* reWin = vtkRenderWindow::New();
    reWin->AddRenderer(renderer);

    vtkRenderWindowInteractor* iren = vtkRenderWindowInteractor::New();
    iren->SetRenderWindow(reWin);
    vtkInteractorStyleTrackballActor* style = vtkInteractorStyleTrackballActor::New();
    iren->SetInteractorStyle(style);

    renderer->SetActiveCamera(camera);
    renderer->ResetCamera();
    renderer->SetBackground(1, 1, 1);

    reWin->SetSize(300, 300);

    reWin->Render();
    iren->Initialize();
    iren->Start();

    points->Delete();
    v_cube->Delete();
    cubeMapper->Delete();
    renderer->Delete();
    reWin->Delete();
    iren->Delete();
    polys->Delete();
#endif
    //opt.close();
    if (segCompleted.size() == 0)segCompleted.emplace_back(0U);

    if (print)
    {
        //這不應該在這裡吧，會害後面分不同 Sensor 顏色輸出時壞掉...
        //pcl::PassThrough<pcl::PointXYZ> pass;
        //pass.setInputCloud(cloud.makeShared());
        //pass.setFilterFieldName("z");
        //pass.setFilterLimits(0.0, 50.0);
        //pass.filter(cloud);
        //pcl::io::savePCDFileBinaryCompressed<pcl::PointXYZ>("simulated_cylinder.pcd", cloud);

        float sAverage = (float)sPts / (float)sCube;
        float sStandard = sqrt((float)sVar / (float)sCube - sAverage * sAverage);
        std::cout << "------------------------------------------------------------------------------\nStatistics\n";
        std::cout << "\nSum: " << sPts << "(" << std::fixed << setprecision(2) << (float)sPts / cloud_size * 100.0f
            << "%), " << sPts / (r * (2 + EIGEN_PI) * seg * 10000.0f) << "pts/cm^2\n"
            << "Cubes: " << sCube << "(Zero: " << sCube - sN << "(" << (float)(sCube - sN) / (float)sCube * 100.0f << "%), N: "
            << sN << "(" << (float)sN / (float)sCube * 100.0f << "%))";
        std::cout << "\n£g = " << std::fixed << setprecision(3) << sAverage << ", £m = " << sStandard;
        std::cout << "\nnCompleted = " << std::fixed << setprecision(1) << nCompleted * precision << endl;
        for (auto i = segCompleted.begin(); i < segCompleted.end();)
        {
            auto it = check_continue(i, segCompleted);
            std::cout << setw(4) << *i * precision << " - " << setw(4) << *it * precision << "\n";
            i = it + 1;
        }
        //std::cout << "------------------------------------------------------------------------------\n";
        //std::cout << setw(6) << "z" << setw(10) << "ave" << setw(10) << "var" << setw(10) << "max" << setw(10) << "min" << endl;
        //for (size_t i = 0; i < vec_max_nPts.size(); i++)
        //{
        //  std::cout << setw(5) << (float)i * precision << "m" << setw(10) << vec_ring_ave.at(i) << setw(10) << vec_ring_var.at(i) << setw(10) << vec_max_nPts.at(i) << setw(10) << vec_min_nPts.at(i) << "\n";
        //}
        //std::cout << "------------------------------------------------------------------------------\n\n";

        //std::fstream opt("diff_amount.csv", std::ios::app | std::ios::out);
        //if (!opt.good())std::cerr << "file failed!\n";
        //opt << "\n, Sum, " << sPts << ", " << (float)sPts / cloud_size * 100.0f << "%, "
        //  << sPts / (r * (2.0f + EIGEN_PI) * seg * 10000.0f) << "pts/cm^2, \n"
        //  << ", Cubes, " << sCube << ", Zero, " << sCube - sN << ", " << (float)(sCube - sN) / sCube * 100.0f << "%, N, "
        //  << sN << ", " << (float)sN / sCube * 100.0f << "%, \n"
        //  << ", £g, " << sAverage << ", £m, " << sStandard << ", \n"
        //  << ", nCompleted, " << nCompleted * precision << *segCompleted.cbegin() * precision << ", \n";
        //opt.close();
    }
    delete[] cube;
    cout << "\n# points: " << cloud.size();
    return (float)nCompleted * precision;
}

Horseshoe::Horseshoe() {}

unsigned Horseshoe::horseshoe(std::vector<Eigen::Vector3d>& v, unsigned length) {
    cloud.width = length - 1 + cloud.size();
    cloud.height = 1;
    //cloud.reserve(cloud.width * cloud.height);

    float ox = v[0].x(), oy = v[0].y(), oz = v[0].z();
    unsigned intersection = 0;
#pragma omp declare reduction (merge : pcl::PointCloud<pcl::PointXYZ> : omp_out.insert(omp_out.end(), omp_in.begin(), omp_in.end()))
#pragma omp parallel for reduction(merge: cloud) reduction(+:intersection) schedule(dynamic, length/THREADS)
    for (unsigned i = 1; i < length; i++)
    {
        float tmpX = 0.0f, tmpY = 0.0f, tmpZ = 0.0f;
        float bx = v[i].x(), by = v[i].y(), bz = v[i].z();
        //At^2+Bt+C, x^2+(y-h)^2=r^2	半圓 
        float A = (bx - ox) * (bx - ox) + (by - oy) * (by - oy);
        float B = 2 * (ox * (bx - ox) + (oy - h) * (by - oy));
        float C = ox * ox + (oy - h) * (oy - h) - r * r;
        float D = B * B - 4 * A * C;
        if (D >= 0.0f)
        {
            float t1R = (-B + sqrt(D)) / (A * 2)/*, t2 = (-BR - sqrt(DR)) / (AR * 2)*/;
            float tmpX1 = ox + (bx - ox) * t1R/*, tmpX2 = ox + (bx - ox) * t2*/;
            float tmpY1 = oy + (by - oy) * t1R/*, tmpY2 = oy + (by - oy) * t2*/;
            float tmpZ1 = oz + (bz - oz) * t1R/*, tmpZ2 = oz + (bz - oz) * t2*/;
            if (tmpY1 >= h)
            {
                tmpX = tmpX1;
                tmpY = tmpY1;
                tmpZ = tmpZ1;            
                auto a = sqrt(pow(tmpX - ox, 2) + pow(tmpY - oy, 2) + pow(tmpZ - oz, 2));//o到b距離
                auto b = sqrt(pow(tmpX - 0, 2) + pow(tmpY - (l - r), 2));//法向量 點到半圓圓心
                auto c = abs((tmpX - ox) * (tmpX - 0) + (tmpY - oy) * (tmpY - (l - r)));
                auto dis = -658.45 * acos(c / (a * b)) + 1041.8;
                if (dis >= a)//判斷角度 y = -658.45x + 1041.8 (x=弧度)  可測最長距離 > o到b距離
                {
                    pcl::PointXYZ tmp{ tmpX,tmpY,tmpZ };
                    cloud.emplace_back(tmp);
                    intersection++;
                }
            }
        }
        float D1 = by - oy;
        if (D1 != 0)
        {
            float t1 = 0 - oy / D1;//y=0
            float tmpX1 = ox + (bx - ox) * t1;
            float tmpY1 = oy + (by - oy) * t1;
            float tmpZ1 = oz + (bz - oz) * t1;
            if (tmpX1 >= -r && tmpX1 <= r)//y=0
            {
                tmpX = tmpX1;
                tmpY = 0;
                tmpZ = tmpZ1;
                auto a = sqrt(pow(tmpX - ox, 2) + pow(tmpY - oy, 2) + pow(tmpZ - oz, 2));//o到b距離
                auto b = 1;//法向量長度
                auto c = abs((tmpY - oy) * 1);
                auto dis = -658.45 * acos(c / (a * b)) + 1041.8;
                if (dis >= a)
                {
                    pcl::PointXYZ tmp{ tmpX,tmpY,tmpZ };
                    cloud.emplace_back(tmp);
                    intersection++; 
                }
            }
        }
        D1 = bx - ox;
        if (D1 != 0) {
            float t1 = (r - ox) / D1;//x=r
            float tmpX1 = ox + (bx - ox) * t1;
            float tmpY1 = oy + (by - oy) * t1;
            float tmpZ1 = oz + (bz - oz) * t1;
            if (tmpY1 >= 0 && tmpY1 <= h) {
                tmpX = r;
                tmpY = tmpY1;
                tmpZ = tmpZ1;
                auto a = sqrt(pow(tmpX - ox, 2) + pow(tmpY - oy, 2) + pow(tmpZ - oz, 2));//o到b距離
                auto b = 1;//法向量長度
                auto c = abs((tmpX - ox) * 1);
                auto dis = -658.45 * acos(c / (a * b)) + 1041.8;
                if (dis >= a)
                {
                    pcl::PointXYZ tmp{ tmpX,tmpY,tmpZ };
                    cloud.emplace_back(tmp);
                    intersection++;
                }             
            }
            float t2 = (-r - ox) / D1;//x=-r
            float tmpX2 = ox + (bx - ox) * t2;
            float tmpY2 = oy + (by - oy) * t2;
            float tmpZ2 = oz + (bz - oz) * t2;
            if (tmpY2 >= 0 && tmpY2 <= h) {
                tmpX = -r;
                tmpY = tmpY2;
                tmpZ = tmpZ2;
                auto a = sqrt(pow(tmpX - ox, 2) + pow(tmpY - oy, 2) + pow(tmpZ - oz, 2));//o到b距離
                auto b = 1;//法向量長度
                auto c = abs((tmpY - oy) * 1);
                auto dis = -658.45 * acos(c / (a * b)) + 1041.8;
                if (dis >= a)
                {
                    pcl::PointXYZ tmp{ tmpX,tmpY,tmpZ };
                    cloud.emplace_back(tmp);
                    intersection++;
                }
            }
        }
    }
    // cout << "\n" << __func__ << ":" << intersection << ", cloud-size: " << cloud.size();
    return intersection;
}

float Horseshoe::statistics(bool print)//統計
{
    cloud_cube.clear();
    Cross_Section::cube_x = (unsigned)(r * 2.0f * scale) + 1U;
    Cross_Section::cube_y = (unsigned)(l * scale) + 1U;
    Cross_Section::cube_z = (unsigned)(seg * scale) + 1U;
    Cross_Section::nCube = cube_x * cube_y * cube_z;
    cloud_cube.reserve(nCube);
    Cross_Section::cube = new unsigned[nCube];
    //cout << "cube_x=" << cube_x << "   cube_y=" << cube_y << "   cube_z=" << cube_z;
    unsigned cloud_size = cloud.size();
    if (print)
    {
        std::cout << "------------------------------------------------------------------------------\n"
            << "Number of points = " << cloud_size << endl;
        std::cout << "Number of cubes = " << cube_x << " * " << cube_y << " * " << cube_z << " = " << nCube
            << "\nSegmentation: " << std::fixed << setprecision(2) << end - start << "m, precision: " << precision << "m\n";
    }

    #pragma omp parallel for num_threads(THREADS) schedule(dynamic, nCube/THREADS)
    for (unsigned i = 0; i < nCube; i++)
    {
        cube[i] = 0U;
    }

    auto start_3 = std::chrono::steady_clock::now();
    #pragma omp parallel for num_threads(THREADS) schedule(dynamic, cloud_size/THREADS)
    for (auto i = 0; i < cloud_size; i++)
    {
        if (cloud.points.at(i).z <= end && cloud.points.at(i).z >= start)
        {
            unsigned _x = (unsigned)((cloud.points.at(i).x + r) * scale);//方塊陣列索引值
            unsigned _y = (unsigned)(cloud.points.at(i).y * scale);
            unsigned _z = (unsigned)((cloud.points.at(i).z - start) * scale);
            //cout << "ncube=" << nCube << " _z * cube_x * cube_y + _x * cube_y + _y =" << _z * cube_x * cube_y + _x * cube_y + _y<<endl;
            ++cube[_z * cube_x * cube_y + _x * cube_y + _y];
        }
    }
    auto end_3 = std::chrono::duration<float>(std::chrono::steady_clock::now() - start_3);
    t_p += end_3.count();

    auto sCube = 0, sPts = 0, sN = 0, sVar = 0;

    float threshold = precision * sqrt(2.0) * 0.9f;//臨界點
    unsigned startScale = start * scale, endScale = end * scale;
    unsigned nCompleted = 0, countCompleted = 0;
    std::vector<unsigned> segCompleted, vec_max_nPts, vec_min_nPts;
    std::vector<float> vec_ring_ave, vec_ring_var;
    vec_max_nPts.reserve(cube_z);
    vec_min_nPts.reserve(cube_z);
    vec_ring_ave.reserve(cube_z);
    vec_ring_var.reserve(cube_z);
    segCompleted.reserve(cube_z);
    //std::fstream opt("var.csv", std::ios::app | std::ios::out);
    //vtkRenderer* renderer = vtkRenderer::New();
    //float corner[8][3] = { { 0, 0, 0 }, { precision, 0, 0 }, { precision, 0, precision }, { 0, 0, precision }, { 0, precision, 0 }, { precision, precision, 0 }, { precision, precision, precision }, { 0, precision, precision } };
    //vtkIdType pts[6][4] = { { 0, 1, 2, 3 }, { 4, 5, 6, 7 }, { 2, 3, 7, 6 }, { 1, 0, 4, 5 }, { 0, 4, 7, 3 }, { 1, 2, 6, 5 } };

    for (unsigned k = 0; k < cube_z; k++)
    {
        float z_ = (float)(k + start * scale) * precision;
        unsigned completed = 0;
        unsigned cross_section_cube = 0;
        unsigned max_nPts = 0;
        unsigned min_nPts = UINT_MAX;
        float ringAverage = 0.0f, ringVar = 0.0f;
        for (unsigned i = 0; i < cube_x; i++)
        {
            float x_ = (float)(i - r * scale) * precision;//變回-r到r
            for (unsigned j = 0; j < cube_y; j++)
            {
                float y_ = (float)j * precision;
                float r_after = (x_ < 0.0f)
                    ? sqrt((x_ + precision) * (x_ + precision) + (y_ - h) * (y_ - h))
                    : sqrt(x_ * x_ + (y_ - h) * (y_ - h));
                if ((int&)y_ >= (int&)h && ((int&)r >= (int&)r_after && abs(r - r_after) - threshold <= eps)) {
                    if (print) {
                        pcl::PointXYZ tmp{ x_,y_,z_ };
                        cloud_cube.emplace_back(tmp);
                    }
                    ++sCube;
                    ++cross_section_cube;
                    auto n = (unsigned long long)cube[k * cube_x * cube_y + i * cube_y + j];
                    if (n > max_nPts)max_nPts = n;
                    if (n < min_nPts)min_nPts = n;
                    if (n)
                    {
                        sVar += n * n;
                        sPts += n;
                        ++sN;
                        ++completed;
                    }
                }
                else if ((int&)y_ <= eps)
                {
                    if (print) {
                        pcl::PointXYZ tmp{ x_,y_,z_ };
                        cloud_cube.emplace_back(tmp);
                    }
                    ++sCube;
                    ++cross_section_cube;
                    auto n = (unsigned long long)cube[k * cube_x * cube_y + i * cube_y + j];
                    //opt << n << ", ";
                    if (n)
                    {
                        sVar += n * n;
                        sPts += n;
                        ++sN;
                        ++completed;
                    }
                    //continue;
                }
                else if (abs(x_) >= r - eps && y_ < h)
                {
                    if (print) {
                        pcl::PointXYZ tmp{ x_,y_,z_ };
                        cloud_cube.emplace_back(tmp);
                    }
                    ++sCube;
                    ++cross_section_cube;
                    auto n = (unsigned long long)cube[k * cube_x * cube_y + i * cube_y + j];
                    if (n)
                    {
                        sVar += n * n;
                        sPts += n;
                        ++sN;
                        ++completed;
                    }
                    //continue;
                }
            }
        }
        //opt << "\n";
        ringAverage /= cross_section_cube;
        float ringVariance = sqrt(ringVar / (float)cross_section_cube - ringAverage * ringAverage);
        vec_ring_ave.emplace_back(ringAverage);
        vec_ring_var.emplace_back(ringVariance);
        vec_max_nPts.emplace_back(max_nPts);
        vec_min_nPts.emplace_back(min_nPts);
        cout << (float)completed / (float)cross_section_cube * 100 <<"%" << endl;
        if ((float)completed / (float)cross_section_cube >= 0.9f)
        {
            ++countCompleted;
            segCompleted.emplace_back((unsigned)(z_ * scale));
        }
        else
        {
            countCompleted = 0U;
        }
        if (countCompleted > nCompleted)nCompleted = countCompleted;
    }
#ifdef VTK
    vtkCamera* camera = vtkCamera::New();
    camera->SetPosition(0, 0, -10);
    camera->SetFocalPoint(0, 0, 0);


    vtkRenderWindow* reWin = vtkRenderWindow::New();
    reWin->AddRenderer(renderer);

    vtkRenderWindowInteractor* iren = vtkRenderWindowInteractor::New();
    iren->SetRenderWindow(reWin);
    vtkInteractorStyleTrackballActor* style = vtkInteractorStyleTrackballActor::New();
    iren->SetInteractorStyle(style);

    renderer->SetActiveCamera(camera);
    renderer->ResetCamera();
    renderer->SetBackground(1, 1, 1);

    reWin->SetSize(300, 300);

    reWin->Render();
    iren->Initialize();
    iren->Start();

    points->Delete();
    v_cube->Delete();
    cubeMapper->Delete();
    renderer->Delete();
    reWin->Delete();
    iren->Delete();
    polys->Delete();
#endif
    //opt.close();
    if (segCompleted.size() == 0)segCompleted.emplace_back(0U);

    if (print)
    {
        float sAverage = (float)sPts / (float)sCube;
        float sStandard = sqrt((float)sVar / (float)sCube - sAverage * sAverage);
        std::cout << "------------------------------------------------------------------------------\nStatistics\n";
        std::cout << "\nSum: " << sPts << "(" << std::fixed << setprecision(2) << (float)sPts / cloud_size * 100.0f
            << "%), " << sPts / (r * (2 + EIGEN_PI) * seg * 10000.0f) << "pts/cm^2\n"
            << "Cubes: " << sCube << "(Zero: " << sCube - sN << "(" << (float)(sCube - sN) / (float)sCube * 100.0f << "%), N: "
            << sN << "(" << (float)sN / (float)sCube * 100.0f << "%))";
        std::cout << "\n£g = " << std::fixed << setprecision(3) << sAverage << ", £m = " << sStandard;
        std::cout << "\nnCompleted = " << std::fixed << setprecision(1) << nCompleted * precision << endl;
        for (auto i = segCompleted.begin(); i < segCompleted.end();)
        {
            auto it = check_continue(i, segCompleted);
            std::cout << setw(4) << *i * precision << " - " << setw(4) << *it * precision << "\n";
            i = it + 1;
        }
        //std::cout << "------------------------------------------------------------------------------\n";
        //std::cout << setw(6) << "z" << setw(10) << "ave" << setw(10) << "var" << setw(10) << "max" << setw(10) << "min" << endl;
        //for (size_t i = 0; i < vec_max_nPts.size(); i++)
        //{
        //  std::cout << setw(5) << (float)i * precision << "m" << setw(10) << vec_ring_ave.at(i) << setw(10) << vec_ring_var.at(i) << setw(10) << vec_max_nPts.at(i) << setw(10) << vec_min_nPts.at(i) << "\n";
        //}
        //std::cout << "------------------------------------------------------------------------------\n\n";

        //std::fstream opt("diff_amount.csv", std::ios::app | std::ios::out);
        //if (!opt.good())std::cerr << "file failed!\n";
        //opt << "\n, Sum, " << sPts << ", " << (float)sPts / cloud_size * 100.0f << "%, "
        //  << sPts / (r * (2.0f + EIGEN_PI) * seg * 10000.0f) << "pts/cm^2, \n"
        //  << ", Cubes, " << sCube << ", Zero, " << sCube - sN << ", " << (float)(sCube - sN) / sCube * 100.0f << "%, N, "
        //  << sN << ", " << (float)sN / sCube * 100.0f << "%, \n"
        //  << ", £g, " << sAverage << ", £m, " << sStandard << ", \n"
        //  << ", nCompleted, " << nCompleted * precision << *segCompleted.cbegin() * precision << ", \n";
        //opt.close();
    }

    delete[] cube;
    cout << "\n# points: " << cloud.size();
    return (float)nCompleted * precision;
}

Arch::Arch() {}

float Arch::a = Arch::l - Arch::r;
float Arch::b = Arch::a * tan(Arch::alpha);
float Arch::R = sqrt(Arch::a * Arch::a + Arch::b * Arch::b) + Arch::r;
float Arch::ix = Arch::a + Arch::r * cos(Arch::alpha);
float Arch::iy = Arch::r * sin(Arch::alpha);
float Arch::f = Arch::R - Arch::b;

unsigned Arch::arch(std::vector<Eigen::Vector3d>& v, unsigned length) {
    cloud.width = length - 1 + cloud.size();
    cloud.height = 1;
    cloud.reserve(cloud.width * cloud.height);

    float ox = v[0].x(), oy = v[0].y(), oz = v[0].z();
    unsigned intersection = 0;
    for (unsigned i = 1; i <= length; i++)
    {
        float tmpX = 0.0f, tmpY = 0.0f, tmpZ = 0.0f;
        float bx = v[i].x(), by = v[i].y(), bz = v[i].z();
        //¥ýºâ¸ò¤j¶ê¥æÂI //??
        auto aR = 0.0f;
        auto bR = -b;
        //Ax+By+C, (x-a)^2+(y-b)^2=r^2
        float AR = (bx - ox) * (bx - ox) + (by - oy) * (by - oy);
        float BR = 2 * ((ox - aR) * (bx - ox) + (oy - bR) * (by - oy));
        float CR = (ox - aR) * (ox - aR) + (oy - bR) * (oy - bR) - R * R;
        float DR = BR * BR - 4 * AR * CR;

        if (DR >= 0.0f)
        {
            float t1R = (-BR + sqrt(DR)) / (AR * 2)/*, t2 = (-BR - sqrt(DR)) / (AR * 2)*/;
            float tmpX1R = ox + (bx - ox) * t1R/*, tmpX2 = ox + (bx - ox) * t2*/;
            float tmpY1R = oy + (by - oy) * t1R/*, tmpY2 = oy + (by - oy) * t2*/;
            float tmpZ1R = oz + (bz - oz) * t1R/*, tmpZ2 = oz + (bz - oz) * t2*/;
            if (tmpY1R > 0.0f /*&& tmpZ1 >= 0*/)
            {
                tmpX = tmpX1R;
                tmpY = tmpY1R;
                tmpZ = tmpZ1R;
            }
            else if (tmpY1R < eps)
            {
                float t = oy / (oy - by);
                tmpX = ox + (bx - ox) * t;
                tmpY = 0.0f;
                tmpZ = oz + (bz - oz) * t;
            }

            if (tmpY < iy /*&& abs(tmpX) > ix*/)
            {
                auto ar = (tmpX > 0) ? a : -a;
                auto br = 0;
                float Ar = (bx - ox) * (bx - ox) + (by - oy) * (by - oy);
                float Br = 2 * ((ox - ar) * (bx - ox) + (oy - br) * (by - oy));
                float Cr = (ox - ar) * (ox - ar) + (oy - br) * (oy - br) - r * r;
                float Dr = Br * Br - 4 * Ar * Cr;
                if (Dr > 0)
                {
                    float t1r = (-Br + sqrt(Dr)) / (Ar * 2)/*, t2 = (-Br - sqrt(Dr)) / (Ar * 2)*/;
                    float tmpX1r = ox + (bx - ox) * t1r/*, tmpX2 = ox + (bx - ox) * t2*/;
                    float tmpY1r = oy + (by - oy) * t1r/*, tmpY2 = oy + (by - oy) * t2*/;
                    float tmpZ1r = oz + (bz - oz) * t1r/*, tmpZ2 = oz + (bz - oz) * t2*/;
                    if (tmpY1r > 0)
                    {
                        tmpX = tmpX1r;
                        tmpY = tmpY1r;
                        tmpZ = tmpZ1r;
                    }
                    else if (tmpY1r < eps)
                    {
                        float t = oy / (oy - by);
                        tmpX = ox + (bx - ox) * t;
                        tmpY = 0.0f;
                        tmpZ = oz + (bz - oz) * t;
                    }
                }
            }
            pcl::PointXYZ tmp{ tmpX,tmpY,tmpZ };
            cloud.emplace_back(tmp);
            ++intersection;
        }
    }
    // cout << "\n" << __func__ << ":" << intersection << ", cloud-size: " << cloud.size();
    return intersection;
}

float Arch::statistics(bool print)
{
    cloud_cube.clear();
    cloud_cube.reserve((cube_x + cube_y) * 2 * cube_z);
    Cross_Section::cube_x = (unsigned)(l * 2.0f * scale) + 1U;
    Cross_Section::cube_y = (unsigned)(f * scale) + 1U;
    Cross_Section::cube_z = (unsigned)(seg * scale) + 1U;
    Cross_Section::nCube = cube_x * cube_y * cube_z;
    Cross_Section::cube = new unsigned[nCube];

    unsigned cloud_size = cloud.size();
    if (print)
    {
        std::cout << "------------------------------------------------------------------------------\n"
            << "Number of points = " << cloud_size << endl;
        std::cout << "Number of cubes = " << cube_x << " * " << cube_y << " * " << cube_z << " = " << nCube
            << "\nSegmentation: " << std::fixed << setprecision(2) << seg - start << "m, precision: " << precision << "m\n";
    }

    #pragma omp parallel for num_threads(THREADS) schedule(dynamic, nCube/THREADS)
    for (unsigned i = 0; i < nCube; i++)
    {
        cube[i] = 0;
    }

    auto start_2 = std::chrono::steady_clock::now();
    #pragma omp parallel for num_threads(THREADS) schedule(dynamic, cloud_size/THREADS)
    for (size_t i = 0; i < cloud_size; i++)
    {
        //if (cloud.points.at(i).z > end || cloud.points.at(i).z < start)continue;
        if (cloud.points.at(i).z <= end && cloud.points.at(i).z >= start)
        {
            unsigned _x = (unsigned)((cloud.points.at(i).x + l) * scale);
            unsigned _y = (unsigned)(cloud.points.at(i).y * scale);
            unsigned _z = (unsigned)((cloud.points.at(i).z - start) * scale);
            ++cube[_z * cube_x * cube_y + _x * cube_y + _y];
        }
    }
    auto end_2 = std::chrono::duration<float>(std::chrono::steady_clock::now() - start_2);
    t_p += end_2.count();


    auto sCube = 0ULL, sPts = 0ULL, sN = 0ULL, sVar = 0ULL;

    float threshold = precision * sqrt(2.0f);
    unsigned startScale = start * scale, endScale = end * scale;
    unsigned nCompleted = 0U, countCompleted = 0U;
    std::vector<unsigned> segCompleted;
    segCompleted.reserve(cube_z);
    //std::vector<std::vector<std::vector<float>>>corner;
    //corner.reserve(nCube);
    //for (size_t i = 0; i < nCube; i++)
    //{
    //  corner[i].reserve(8);
    //}
    //for (size_t i = 0; i < nCube; i++)
    //{
    //  for (size_t j = 0; j < 8; j++)
    //  {
    //      corner[i][j].reserve(3);
    //  }
    //}
    for (unsigned k = 0; k < cube_z; k++)
    {
        float z_ = (float)(k + start * scale) * precision;
        unsigned completed = 0U;
        unsigned cross_section_cube = 0U;
        for (unsigned i = 0; i < cube_x; i++)
        {
            float x_ = (float)(i - l * scale) * precision;
            for (unsigned j = 0; j < cube_y; j++)
            {
                float y_ = (float)j * precision;
                if (abs(y_) <= eps)
                {
                    if (print) {
                        pcl::PointXYZ tmp{ x_,y_,z_ };
                        cloud_cube.emplace_back(tmp);
                    }
                    ++sCube;
                    ++cross_section_cube;
                    auto index = k * cube_x * cube_y + i * cube_y + j;
                    auto n = (unsigned long long)cube[index];
                    if (n)
                    {
                        //for (size_t px = 0; px < 2; px++)
                        //{
                        //  for (size_t py = 0; py < 2; py++)
                        //  {
                        //      for (size_t pz = 0; pz < 2; pz++)
                        //      {
                        //          auto np = (px + 1) * (py + 1) * (pz + 1) - 1;
                        //          corner[index][np][0] = x_ + px * 0.1;
                        //          corner[index][np][1] = y_ + py * 0.1;
                        //          corner[index][np][2] = z_ + pz * 0.1;
                        //      }
                        //  }
                        //}

                        sVar += n * n;
                        sPts += n;
                        ++sN;
                        ++completed;
                    }
                    continue;
                }

                float r_cube = 0.0;
                float r_after = 0.0;
                bool isEmpty = true;
                if (abs(x_) <= ix)
                {
                    r_after = (x_ < 0)
                        ? sqrt((x_ + precision) * (x_ + precision) + (y_ + b) * (y_ + b))
                        : sqrt(x_ * x_ + (y_ + b) * (y_ + b));
                    if ((int&)R >= (int&)r_after && abs(R - r_after) - threshold <= eps)isEmpty = false;
                }
                else if (abs(x_) <= l)
                {
                    auto aS = (x_ > ix) ? a : -a;
                    r_after = (x_ < 0)
                        ? sqrt((x_ - aS + precision) * (x_ - aS + precision) + y_ * y_)
                        : sqrt((x_ - aS) * (x_ - aS) + y_ * y_);
                    if ((int&)r >= (int&)r_after && abs(r - r_after) - threshold <= eps)isEmpty = false;
                }

                if (!isEmpty)
                {
                    if (print) {
                        pcl::PointXYZ tmp{ x_,y_,z_ };
                        cloud_cube.emplace_back(tmp);
                    }
                    ++sCube;
                    ++cross_section_cube;
                    auto n = (unsigned long long)cube[k * cube_x * cube_y + i * cube_y + j];
                    if (n)
                    {
                        sVar += n * n;
                        sPts += n;
                        ++sN;
                        ++completed;
                    }
                }
            }
        }

        if ((float)completed / (float)cross_section_cube >= 0.9f)
        {
            ++countCompleted;
            segCompleted.emplace_back((unsigned)(z_ * scale));
        }
        else
        {
            countCompleted = 0;
        }
        if (countCompleted > nCompleted)nCompleted = countCompleted;
    }
    if (segCompleted.size() == 0)segCompleted.emplace_back(0);
    if (print)
    {
        float sAverage = (float)sPts / sCube;
        float sStandard = sqrt((float)sVar / sCube - sAverage * sAverage);
        float surface_area = (r * 2 + a * 2 + 2 * EIGEN_PI * r * 2 * alpha * 180.0f / EIGEN_PI / 360.0f + 2 * EIGEN_PI * R * 2 * (90 - alpha * 180.0f / EIGEN_PI) / 360.0f) * seg * 10000.0;
        std::cout << "------------------------------------------------------------------------------\nStatistics\n";
        std::cout << "\nSum: " << sPts << "(" << std::fixed << setprecision(2) << (float)sPts / cloud_size * 100.0f
            << "%), " << sPts / surface_area << "pts/cm^2\n"
            << "Cubes: " << sCube << "(Zero: " << sCube - sN << "(" << (float)(sCube - sN) / sCube * 100.0f << "%), N: "
            << sN << "(" << (float)sN / sCube * 100.0f << "%))";
        std::cout << "\n£g = " << std::fixed << setprecision(3) << sAverage << ", £m = " << sStandard;
        std::cout << "\nnCompleted = " << std::fixed << setprecision(1) << nCompleted * precision << endl;
        for (auto i = segCompleted.begin(); i < segCompleted.end();)
        {
            auto it = check_continue(i, segCompleted);
            std::cout << setw(4) << *i * precision << " - " << setw(4) << *it * precision << "\n";
            i = it + 1;
        }
        std::cout << "------------------------------------------------------------------------------\n\n";
        //std::fstream opt("diff_amount.csv", std::ios::app | std::ios::out);
        //if (!opt.good())std::cerr << "file failed!\n";
        //opt << "\n, Sum, " << sPts << ", " << (float)sPts / cloud_size * 100.0f << "%, "
        //  << sPts / (r * (2 + EIGEN_PI) * seg * 10000.0) << "pts/cm^2, \n"
        //  << ", Cubes, " << sCube << ", Zero, " << sCube - sN << ", " << (float)(sCube - sN) / sCube * 100.0f << "%, N, "
        //  << sN << ", " << (float)sN / sCube * 100.0f << "%, \n"
        //  << ", £g, " << sAverage << ", £m, " << sStandard << ", \n"
        //  << ", nCompleted, " << nCompleted * precision << *segCompleted.cbegin() * precision << ", \n";
        //opt.close();
    }
    delete[] cube;
    return (float)nCompleted * precision;
}
