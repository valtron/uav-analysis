import os, sys
import csv

input_skeleton_pre = """&aircraft_data
   aircraft%cname          = 'UAV'      ! M  name of aircraft
   aircraft%ctype          = 'SymCPS UAV Design' ! M  type of aircraft
   aircraft%num_wings      = 0 ! M number of wings in aircraft
   aircraft%mass          = 1.251158179559168
   aircraft%x_cm          = 0.004262667602479541
   aircraft%y_cm          = 0.004262257946293246
   aircraft%z_cm          = -12.185845568919374
   aircraft%x_fuse          = 0.004262667602479541
   aircraft%y_fuse          = 0.004262257946293246
   aircraft%z_fuse          = -12.185845568919374
   aircraft%X_fuseuu      = 15197.920416918165
   aircraft%Y_fusevv      = 15197.920931262906
   aircraft%Z_fuseww      = 27152.33937824657
   aircraft%Ixx           = 26102.092848709824
   aircraft%Iyy           = 26102.092914663815
   aircraft%Izz           = 48260.9402281981
   aircraft%uc_initial     = 0.4d0, 0.5d0, 0.6d0, 0.7d0 ! inputs for controls
   aircraft%time           = 0.d0        ! initial time (default = 0.)
   aircraft%dt             = 1.d-03      ! s  fixed time step
   aircraft%dt_output      = 1.0d0       ! s  time between output lines
   aircraft%time_end       = 1000.d0        ! s  end time 
   aircraft%Unwind         = 0.d0        ! North wind speed in world frame
   aircraft%Vewind         = 0.d0        ! East wind speed in  world frame
   aircraft%Wdwind         = 0.d0        ! Down wind speed in world frame
   aircraft%debug          = 0           ! verbose printouts from fderiv
   aircraft%num_propellers  = 1
   aircraft%num_batteries   = 1
   aircraft%i_analysis_type = 1 
   aircraft%x_initial      = 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 1.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0
   aircraft%uc_initial     = 0.5d0, 0.5d0, 0.5d0, 0.5d0

!   Propeller(1 uses components named Prop_2, Motor_2, ESC_2"""

input_skeleton_post = """!   Controls
   control%i_flight_path = 4
   control%requested_lateral_speed = 37.0 
   control%requested_vertical_speed = 0.0
   control%iaileron = 5 
   control%iflap = 6 
   control%Q_position = 1.0 
   control%Q_velocity = 1.0
   control%Q_angular_velocity = 1.0 
   control%Q_angles = 1.0 
   control%R= 50.0 
/
"""


def read_BEMPs(dir):
    with open(os.path.join(dir, 'Motor.csv'), newline='') as csvfile:
        motors = list(csv.DictReader(csvfile))
    with open(os.path.join(dir, 'Propeller.csv'), newline='') as csvfile:
        props = list(csv.DictReader(csvfile))
    with open(os.path.join(dir, 'Battery.csv'), newline='') as csvfile:
        batteries = list(csv.DictReader(csvfile))

    return batteries, motors, props


def generate_input(b, m, p):
    str_return = ""
    str_return += input_skeleton_pre
    str_return += "   propeller(1)%cname = '{}'\n".format(p["Name"])
    str_return += "   propeller(1)%ctype = 'MR'\n"
    str_return += "   propeller(1)%prop_fname = '../flight-dynamics-model-master/propeller/{}'\n".format(
        p["Performance_File"])
    str_return += "   propeller(1)%x = 0\n"
    str_return += "   propeller(1)%y = 0\n"
    str_return += "   propeller(1)%z = 0\n"
    str_return += "   propeller(1)%nx = -0.0\n"
    str_return += "   propeller(1)%ny = -0.0\n"
    str_return += "   propeller(1)%nz = -1.0\n"
    str_return += "   propeller(1)%radius = {}\n".format(float(p["DIAMETER"]) / 2.0)
    str_return += "   propeller(1)%Ir = {}\n".format(10.0)  # unknown value but not used
    str_return += "   propeller(1)%motor_fname = '{}'\n".format(m["MODEL"])  # not used
    str_return += "   propeller(1)%KV = {}\n".format(m["KV"])
    str_return += "   propeller(1)%KT = {}\n".format(m["KT"])
    str_return += "   propeller(1)%I_max = {}\n".format(m["MAX_CURRENT"])
    str_return += "   propeller(1)%I_idle = {}\n".format(m["IO_IDLE_CURRENT@10V"])
    str_return += "   propeller(1)%maxpower = {}\n".format(m["MAX_POWER"])
    str_return += "   propeller(1)%Rw = {}\n".format(float(m["INTERNAL_RESISTANCE"]) / 1000.0)
    str_return += "   propeller(1)%icontrol = 4\n"
    str_return += "   propeller(1)%ibattery = 1\n"
    str_return += "   propeller(1)%spin = -1\n"
    str_return += "\n"
    str_return += "!   Battery(1) is component named: Battery_0\n"
    str_return += "   battery(1)%num_cells = {}\n".format(int(b["NUMBER_OF_CELLS"][0]))
    str_return += "   battery(1)%voltage = {}\n".format(b["VOLTAGE"])
    str_return += "   battery(1)%capacity = {}\n".format(b["CAPACITY"])
    str_return += "   battery(1)%C_Continuous = {}\n".format(b["CONT_DISCHARGE_RATE"])
    str_return += "   battery(1)%C_Peak = {}\n".format(b["PEAK_DISCHARGE_RATE"])
    str_return += "\n"
    str_return += input_skeleton_post

    return str_return


def main():
    batteries, motors, props = read_BEMPs('../uav_analysis/data/')
    output_header = "Battery,Motor,Propeller," \
                    "MV_omega_rad,MV_omega_RPM,MV_Voltage,MV_Thrust,MV_Torque,MV_Power,MV_Current,MV_Efficiency," \
                    "MP_omega_rad,MP_omega_RPM,MP_Voltage,MP_Thrust,MP_Torque,MP_Power,MP_Current,MP_Efficiency," \
                    "MC_omega_rad,MC_omega_RPM,MC_Voltage,MC_Thrust,MC_Torque,MC_Power,MC_Current,MC_Efficiency," \
                    "MaxPower,MaxCur,PeakCur,ContCur\n"
    res_fname = 'results2.csv'
    res_file = open(res_fname, 'w')
    res_file.write(output_header)
    cnt = 1
    num_of_combinations = len(batteries) * len(motors) * len(props)
    print("Progress: ")
    for b in batteries:
        for m in motors:
            for p in props:
                file_object = open('input.txt', 'w')
                file_object.write(generate_input(b, m, p))
                file_object.close()
                combo_name = "B_" + b["MODEL"] + \
                             "_M_" + m["MODEL"] + \
                             "_P_" + p["Performance_File"] + "_out.txt"
                combo_name = "output.csv"
                cmd = "./../flight-dynamics-model-master/bin/new_fdm < input.txt > {}".format(combo_name)
                os.system(cmd)
                file_object = open(combo_name, 'r')
                res = file_object.read()
                file_object.close()
                res = res.split("\n")
                MV = res[3].split()
                MP = res[4].split()
                MC = res[5].split()

                row = ""
                row += b["MODEL"] + ","
                row += m["MODEL"] + ","
                row += p["Performance_File"] + ","
                row += MV[3] + "," + MV[4] + "," + MV[5] + "," + MV[6] + "," + MV[7] + "," + MV[8] + "," + MV[9] + "," + \
                       MV[10] + ","
                row += MP[3] + "," + MP[4] + "," + MP[5] + "," + MP[6] + "," + MP[7] + "," + MP[8] + "," + MP[9] + "," + \
                       MP[10] + ","
                row += MC[3] + "," + MC[4] + "," + MC[5] + "," + MC[6] + "," + MC[7] + "," + MC[8] + "," + MC[9] + "," + \
                       MC[10] + ","
                row += MC[11] + "," + MC[12] + "," + MC[13] + "," + MC[14] + "\n"
                res_file.write(row)

                sys.stdout.write('\r')
                sys.stdout.write("{:10.2f}%".format(100 * cnt / num_of_combinations))
                sys.stdout.flush()
                cnt += 1
    print()
    print("---------------------------------------")
    print("Results are saved in file: ", res_fname)
    res_file.close()


if __name__ == '__main__':
    main()
