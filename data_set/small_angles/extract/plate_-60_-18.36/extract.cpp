// =============================================================================
// PROJECT CHRONO - http://projectchrono.org
//
// Copyright (c) 2019 projectchrono.org
// All rights reserved.
//
// Use of this source code is governed by a BSD-style license that can be found
// in the LICENSE file at the top level of the distribution and at
// http://projectchrono.org/license-chrono.txt.
//
//This demo simulate a plate intrude into a box area of granular materials with certain attack and intrusion angles

#include <iostream>
#include <vector>
#include <string>
#include "chrono/physics/ChBodyEasy.h"
#include "chrono/core/ChGlobal.h"
#include "chrono_thirdparty/filesystem/path.h"
#include "chrono/physics/ChSystemSMC.h"
#include "chrono/physics/ChBody.h"
#include "chrono/physics/ChForce.h"
#include "chrono/utils/ChUtilsSamplers.h"
#include "chrono/timestepper/ChTimestepper.h"
#include "chrono_granular/api/ChApiGranularChrono.h"
#include "chrono_granular/physics/ChGranular.h"
#include "chrono_granular/physics/ChGranularTriMesh.h"
#include "chrono/utils/ChUtilsCreators.h"
#include "chrono_granular/utils/ChGranularJsonParser.h"

using namespace chrono;
using namespace chrono::granular;

void ShowUsage(std::string name) {
    std::cout << "usage: " + name + " <json_file>" << std::endl;
}

void writeMeshFrames(std::ostringstream& outstream, ChBody& body, std::string obj_name, float mesh_scaling) {
    outstream << obj_name << ",";

    // Get frame position
    ChFrame<> body_frame = body.GetFrame_REF_to_abs();
    ChQuaternion<> rot = body_frame.GetRot();
    ChVector<> pos = body_frame.GetPos();

    // Get basis vectors
    ChVector<> vx = rot.GetXaxis();
    ChVector<> vy = rot.GetYaxis();
    ChVector<> vz = rot.GetZaxis();

    // Output in order
    outstream << pos.x() << ",";
    outstream << pos.y() << ",";
    outstream << pos.z() << ",";
    outstream << vx.x() << ",";
    outstream << vx.y() << ",";
    outstream << vx.z() << ",";
    outstream << vy.x() << ",";
    outstream << vy.y() << ",";
    outstream << vy.z() << ",";
    outstream << vz.x() << ",";
    outstream << vz.y() << ",";
    outstream << vz.z() << ",";
    outstream << mesh_scaling << "," << mesh_scaling << "," << mesh_scaling;
    outstream << "\n";
}
const double time_settle = 0.5;
const double time_intrude = 1.0;
const double time_pre = time_settle + time_intrude;
constexpr float F_CGS_TO_SI = 1e-5;
int main(int argc, char* argv[]) {

    std::ofstream out_as("output_plate_forces.csv");
    std::ofstream out_pos("output_plate_positions.csv");
    sim_param_holder params;
    if (argc != 2 || ParseJSON(argv[1], params) == false) {
        ShowUsage(argv[0]);
        return 1;
    }

    float iteration_step = params.step_size;

    ChGranularChronoTriMeshAPI apiSMC_TriMesh(params.sphere_radius, params.sphere_density,
                                              make_float3(params.box_X, params.box_Y, params.box_Z));

    ChSystemGranularSMC_trimesh& gran_sys = apiSMC_TriMesh.getGranSystemSMC_TriMesh();
    double fill_bottom = -params.box_Z / 2.0; // -200/2
    double fill_top = params.box_Z / 2.0;     // 200/4

   // chrono::utils::PDSampler<float> sampler(2.4f * params.sphere_radius);
    chrono::utils::HCPSampler<float> sampler(2.05 * params.sphere_radius);

    // leave a 4cm margin at edges of sampling
    ChVector<> hdims(params.box_X / 2-0.5 , params.box_Y / 2-0.5 , 0);
    ChVector<> center(0, 0, fill_bottom + 2.0 * params.sphere_radius);
    std::vector<ChVector<float>> body_points;

    // Shift up for bottom of box
    center.z() += 3 * params.sphere_radius;
    while (center.z() < fill_top) {
        std::cout << "Create layer at " << center.z() << std::endl;
        auto points = sampler.SampleBox(center, hdims);
        body_points.insert(body_points.end(), points.begin(), points.end());
        center.z() += 2.05 * params.sphere_radius;
    }

    apiSMC_TriMesh.setElemsPositions(body_points);

    gran_sys.set_BD_Fixed(true);
    std::function<double3(float)> pos_func_wave = [&params](float t) {
        double3 pos = {0, 0, 0};

        double t0 = 0.5;
        double freq = CH_C_PI / 4;

        if (t > t0) {
            pos.x = 0.1 * params.box_X * std::sin((t - t0) * freq);
        }
        return pos;
    };

    // gran_sys.setBDWallsMotionFunction(pos_func_wave);

    gran_sys.set_K_n_SPH2SPH(params.normalStiffS2S);
    gran_sys.set_K_n_SPH2WALL(params.normalStiffS2W);
    gran_sys.set_K_n_SPH2MESH(params.normalStiffS2M);

    gran_sys.set_Gamma_n_SPH2SPH(params.normalDampS2S);
    gran_sys.set_Gamma_n_SPH2WALL(params.normalDampS2W);
    gran_sys.set_Gamma_n_SPH2MESH(params.normalDampS2M);

    gran_sys.set_K_t_SPH2SPH(params.tangentStiffS2S);
    gran_sys.set_K_t_SPH2WALL(params.tangentStiffS2W);
    gran_sys.set_K_t_SPH2MESH(params.tangentStiffS2M);

    gran_sys.set_Gamma_t_SPH2SPH(params.tangentDampS2S);
    gran_sys.set_Gamma_t_SPH2WALL(params.tangentDampS2W);
    gran_sys.set_Gamma_t_SPH2MESH(params.tangentDampS2M);

    gran_sys.set_Cohesion_ratio(params.cohesion_ratio);
    gran_sys.set_Adhesion_ratio_S2M(params.adhesion_ratio_s2m);
    gran_sys.set_Adhesion_ratio_S2W(params.adhesion_ratio_s2w);
    gran_sys.set_gravitational_acceleration(params.grav_X, params.grav_Y, params.grav_Z);

    gran_sys.set_fixed_stepSize(params.step_size);
    gran_sys.set_friction_mode(GRAN_FRICTION_MODE::MULTI_STEP);
    gran_sys.set_timeIntegrator(GRAN_TIME_INTEGRATOR::CENTERED_DIFFERENCE);
    gran_sys.set_static_friction_coeff_SPH2SPH(params.static_friction_coeffS2S);
    gran_sys.set_static_friction_coeff_SPH2WALL(params.static_friction_coeffS2W);
    gran_sys.set_static_friction_coeff_SPH2MESH(params.static_friction_coeffS2M);

    std::string mesh_filename(GetChronoDataFile("granular/demo_GRAN_plate/plate.obj"));
    std::vector<string> mesh_filenames(1, mesh_filename);

    std::vector<float3> mesh_translations(1, make_float3(0.f, 0.f, 0.f));

    float ball_radius = 20.f;
    float length = 3.81;
    float width = 2.54;
    float thickness = 0.64;
    std::vector<ChMatrix33<float>> mesh_rotscales(1, ChMatrix33<float>(0.50));

    float plate_density = 2.7;//params.sphere_density / 100.f;
    float plate_mass = (float)length * width * thickness * plate_density ;
    std::vector<float> mesh_masses(1, plate_mass);

    apiSMC_TriMesh.load_meshes(mesh_filenames, mesh_rotscales, mesh_translations, mesh_masses);

    gran_sys.setOutputMode(params.write_mode);
    gran_sys.setVerbose(params.verbose);
    filesystem::create_directory(filesystem::path(params.output_dir));

    unsigned int nSoupFamilies = gran_sys.getNumTriangleFamilies();
    std::cout << nSoupFamilies << " soup families" << std::endl;
    double* meshPosRot = new double[7 * nSoupFamilies];
    float* meshVel = new float[6 * nSoupFamilies]();

    gran_sys.initialize();

// create a plate for simulation
    ChSystemSMC sys_plate;
    sys_plate.SetContactForceModel(ChSystemSMC::ContactForceModel::Hooke);
    sys_plate.SetTimestepperType(ChTimestepper::Type::EULER_EXPLICIT);
    sys_plate.Set_G_acc(ChVector<>(0, 0, -980));
  //  auto rigid_plate = std::make_shared<ChBodyEasyBox>(length, width, thickness, plate_density, true, true);
    std::shared_ptr<ChBody> rigid_plate(sys_plate.NewBody());
    rigid_plate->SetMass(plate_mass);
    rigid_plate->SetPos(ChVector<>(0,0,12));
    float inertiax = 1 / 12 * plate_mass*(thickness* thickness +width*width);
    float inertiay = 1 / 12 * plate_mass * (thickness * thickness + length * length);
    float inertiaz = 1 / 12 * plate_mass * (length * length + width * width);
    rigid_plate->SetInertiaXX(ChVector<>(inertiax, inertiay, inertiaz));
    //sys_plate.AddBody(rigid_plate);
    rigid_plate->SetBodyFixed(true);
    sys_plate.AddBody(rigid_plate); 
    unsigned int out_fps = 50;
    std::cout << "Rendering at " << out_fps << "FPS" << std::endl;

    unsigned int out_steps = (unsigned int)(1.0 / (out_fps * iteration_step));

    int currframe = 0;
    unsigned int curr_step = 0;
    gran_sys.disableMeshCollision();
    clock_t start = std::clock();
    bool plate_released = false;
    bool intrude_state = false;
    double max_z = gran_sys.get_max_z();
    double start_gamma = -0.102*CH_C_PI;
    double resolution = start_gamma / 10;
    double v_intrude = -15.0;
    
    for (int i = 0; i < 1; i++) {
        // define the velocity, body orientations for the plate
        double gamma = start_gamma-resolution*i;
        double belta =-CH_C_PI / 3;
        std::cout<< "gamma="  << gamma * 180 / CH_C_PI << "belta=" << belta * 180 / CH_C_PI << std::endl;
        out_as << gamma*180/CH_C_PI << ',' << belta * 180 / CH_C_PI << '\n';
        out_pos << gamma * 180 / CH_C_PI << ',' <<  belta * 180 / CH_C_PI << '\n';
        double vx = -cos(gamma) * 4;
        double vz = -sin(gamma) * 4;
        
        double start_height = length / 2 * sin(belta);
        int counter = 0;
        rigid_plate->SetPos(ChVector<>(10, 0, 12));
	    rigid_plate->SetBodyFixed(true);
	    plate_released=false;
        intrude_state = false;
        gran_sys.disableMeshCollision();
        for (double t = 0; t < (double)params.time_end; t += iteration_step, curr_step++) {

            if (t >= time_settle && plate_released == false) {
                gran_sys.enableMeshCollision();

                rigid_plate->SetBodyFixed(false);
                max_z = gran_sys.get_max_z();
               // rigid_plate->SetPos(ChVector<>(0, 0, max_z + start_height));
                rigid_plate->SetPos(ChVector<>(10, 0, 12));
                rigid_plate->SetPos_dt(ChVector<>(0, 0, v_intrude));
                rigid_plate->SetRot(Q_from_AngAxis(belta, VECT_Y));
                //   rigid_plate->SetRot(ChQuaternion<>(0.707, 0, 0.707, 0));
                plate_released = true;
                std::cout << "Releasing ball" << std::endl;
            }
            else if (t >= time_settle && t < (time_settle + time_intrude)&&plate_released == true) {
                rigid_plate->SetPos_dt(ChVector<>(0, 0, v_intrude));
                rigid_plate->SetRot(Q_from_AngAxis(belta, VECT_Y));
                //   rigid_plate->SetRot(ChQuaternion<>(0.707, 0, 0.707, 0));
               // std::cout << "Plate intruding" << std::endl;
            }
            else if (t >= (time_settle + time_intrude) && plate_released == true) {
               // intrude_state = true;
                rigid_plate->SetPos_dt(ChVector<>(vx, 0, vz));
                rigid_plate->SetRot(Q_from_AngAxis(belta, VECT_Y));
            }
            auto plate_pos = rigid_plate->GetPos();
            auto plate_rot = rigid_plate->GetRot();
            /*
            auto plate_vel = rigid_plate->GetPos_dt();
            auto plate_ang_vel = rigid_plate->GetWvel_loc();
            plate_ang_vel = rigid_plate->GetRot().GetInverse().Rotate(plate_ang_vel);
            */
            meshPosRot[0] = plate_pos.x();
            meshPosRot[1] = plate_pos.y();
            meshPosRot[2] = plate_pos.z();
            meshPosRot[3] = plate_rot[0];
            meshPosRot[4] = plate_rot[1];
            meshPosRot[5] = plate_rot[2];
            meshPosRot[6] = plate_rot[3];

            meshVel[0] = (float)vx;//plate_vel.x();
            meshVel[1] = (float)0;// plate_vel.y();
            meshVel[2] = (float)vz; //plate_vel.z();
            meshVel[3] = (float)0;// plate_ang_vel.x();
            meshVel[4] = (float)0;// plate_ang_vel.y();
            meshVel[5] = (float)0;// plate_ang_vel.z();

            gran_sys.meshSoup_applyRigidBodyMotion(meshPosRot, meshVel);

            gran_sys.advance_simulation(iteration_step);
            sys_plate.DoStepDynamics(iteration_step);

            float plate_force[6];
            gran_sys.collectGeneralizedForcesOnMeshSoup(plate_force);

            rigid_plate->Empty_forces_accumulators();
            rigid_plate->Accumulate_force(ChVector<>(plate_force[0], plate_force[1], plate_force[2]), plate_pos, false);
            rigid_plate->Accumulate_torque(ChVector<>(plate_force[3], plate_force[4], plate_force[5]), false);
            /*       std::cout <<rigid_plate->GetPos()[2]<<','<< rigid_plate->Get_accumulated_force()[0] * F_CGS_TO_SI << ','
                             << rigid_plate->Get_accumulated_force()[1] * F_CGS_TO_SI << ','
                             << rigid_plate->Get_accumulated_force()[2] * F_CGS_TO_SI <<','<<gran_sys.get_max_z()<<','<<gran_sys.getNumSpheres()<< std::endl;
                   out_as << t << "," << rigid_plate->Get_accumulated_force()[0] * F_CGS_TO_SI << ","
                          << rigid_plate->Get_accumulated_force()[1] * F_CGS_TO_SI<<","
                          << rigid_plate->Get_accumulated_force()[2] * F_CGS_TO_SI
                          << '\n';
                          */
            if (counter % 100 == 0) {
                std::cout << t << ',' << rigid_plate->GetPos()[2] << ',' << plate_force[0] * F_CGS_TO_SI << ','
                    << plate_force[1] * F_CGS_TO_SI << ','
                    << plate_force[2] * F_CGS_TO_SI << ',' << gran_sys.get_max_z() << ',' << gran_sys.getNumSpheres() << std::endl;
                out_as << t << "," << plate_force[0] * F_CGS_TO_SI << ","
                    << plate_force[1] * F_CGS_TO_SI << ","
                    << plate_force[2] * F_CGS_TO_SI
                    << '\n';
                out_pos << t << "," << rigid_plate->GetPos()[0] << "," << rigid_plate->GetPos()[1] << ","
                    << rigid_plate->GetPos()[2] << "," << gran_sys.get_max_z() << "\n";
            }
            if (curr_step % out_steps == 0) {
                std::cout << "Rendering frame " << currframe << std::endl;
                /*     char filename[100];
                     sprintf(filename, "%s/step%06d", params.output_dir.c_str(), currframe++);
                     gran_sys.writeFile(std::string(filename));

                     std::string mesh_output = std::string(filename) + "_meshframes.csv";
                     std::ofstream meshfile(mesh_output);
                     std::ostringstream outstream;
                     outstream << "mesh_name,dx,dy,dz,x1,x2,x3,y1,y2,y3,z1,z2,z3,sx,sy,sz\n";
                     writeMeshFrames(outstream, *ball_body, mesh_filename, ball_radius);
                     meshfile << outstream.str();*/
            }
            counter++;
        }
    }
    clock_t end = std::clock();
    double total_time = ((double)(end - start)) / CLOCKS_PER_SEC;

    std::cout << "Time: " << total_time << " seconds" << std::endl;

    delete[] meshPosRot;
    delete[] meshVel;
    out_as.close();
    out_pos.close();
    return 0;
}
