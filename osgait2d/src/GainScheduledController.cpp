#include <iostream>
#include <vector>
#include <stdexcept>

#include <OpenSim/OpenSim.h>

double ComputePercentGait (const double current_time, const SimTK::Vector& heel_strike_times) {
    // Given a vector of monotonically increasing heel strike times and the
    // current time, return the percent of the gait cycle of the current time.

    int length = heel_strike_times.size();

    // current_time has to be in between the first and last heel strike times
    if (current_time < heel_strike_times(0) || current_time > heel_strike_times(length - 1))
    {
        throw std::invalid_argument("The heel strikes times must bound the current time.");
    }

    int i;
    for (i = 0; i < length; i++){
        if (heel_strike_times(i) > current_time){
            break;
        }
    }

    double time_at_0 = heel_strike_times(i - 1);
    double time_at_100 = heel_strike_times(i);

    return (current_time - time_at_0) / (time_at_100 - time_at_0);
}

double interpolate(double n1, double n3, double d1, double d2, double d3){
    // n3 - n1   n2 - n1           d2 - d1
    // ------- = -------  => n2 =  ------- (n3 - n1) + n1
    // d3 - d1   d2 - d1           d3 - d1

    return (d2 - d1) / (d3 - d1) * (n3 - n1) + n1;

}

double extrapolate(double n1, double n2, double d1, double d2, double d3){
    // n3 - n1   n2 - n1         (d3 - d1)
    // ------- = ------- => n3 = --------- (n2 - n1) + n1
    // d3 - d1   d2 - d1         (d2 - d1)

    return (d3 - d1) / (d2 - d1) * (n2 - n1) + n1;
}

SimTK::Matrix InterpolateGainArray(const double percent_gait_cycle, const std::vector<SimTK::Matrix>& gain_array){
    // Given a n x q x p gain array (n : num percent gain discretizations, q: num
    // controls, p : num sensors) and a percent gait cycle value from 0.0 to 1.0,
    // return a gain matrix which is a linear interpolation of the matrix elements
    // of the two adjacent matrices.


    if (percent_gait_cycle < 0.0 || percent_gait_cycle > 1.0)
    {
        throw std::invalid_argument("The percent gait cycle must be between 0.0 and 1.0.");
    }

    // Make a percent gain vector for this gain array based on n.
    int n = gain_array.size();
    SimTK::Vector percent_gain_vec(n);
    for (int i = 0; i < n; i++){
        percent_gain_vec(i) = i * 1.0 / n;
    }
    std::cout << "percent_gait_vec = " << percent_gain_vec << std::endl;

    // Find the indice which comes after the provided percent gait cycle.
    int i;
    for (i = 0; i < n; i++){
        if (percent_gain_vec(i) > percent_gait_cycle){
            break;
        }
    }

    bool need_to_extrapolate = i == n;

    std::cout << "i = " << i << std::endl;

    double first_percent;
    double second_percent;

    SimTK::Matrix first_gain_matrix;
    SimTK::Matrix second_gain_matrix;

    if (need_to_extrapolate)
    {
        first_percent = percent_gain_vec(i - 2);
        second_percent = percent_gain_vec(i - 1);
        first_gain_matrix = gain_array[i - 2];
        second_gain_matrix = gain_array[i - 1];
    }
    else
    {
        first_percent = percent_gain_vec(i - 1);
        second_percent = percent_gain_vec(i);
        first_gain_matrix = gain_array[i - 1];
        second_gain_matrix = gain_array[i];
    }

    std::cout << "first_percent = " << first_percent << std::endl;
    std::cout << "second_percent = " << second_percent << std::endl;
    std::cout << "known_percent = " << percent_gait_cycle << std::endl;

    std::cout << "first_gain_matrix = " << first_gain_matrix << std::endl;
    std::cout << "second_gain_matrix = " << second_gain_matrix << std::endl;

    int q = first_gain_matrix.nrow();
    int p = first_gain_matrix.ncol();

    SimTK::Matrix interpolated_gain_matrix(q, p);

    for (int i = 0; i < q; i++)
    {
        for (int j = 0; j < p; j++)
        {
            if (need_to_extrapolate)
            {
                interpolated_gain_matrix(i, j) = extrapolate(first_gain_matrix(i, j),
                                                             second_gain_matrix(i, j),
                                                             first_percent,
                                                             second_percent,
                                                             percent_gait_cycle);
            }
            else
            {

                interpolated_gain_matrix(i, j) = interpolate(first_gain_matrix(i, j),
                                                             second_gain_matrix(i, j),
                                                             first_percent,
                                                             percent_gait_cycle,
                                                             second_percent);
            }
         }
     }

    return interpolated_gain_matrix;
}

class GainScheduledController : public OpenSim::Controller {
OpenSim_DECLARE_CONCRETE_OBJECT(GainScheduledController, OpenSim::Controller);

public:
    // This controller needs three pieces of information on construction:
    // gain_array : n vector of q x p matrices where each matrix is gain matrix
    // corresponding to a percentage of the gait cycle
    // t_star_array : n vector of q x 1 matrices where each matrix is T* vector
    // corresponding to a percentage of the gait cycle
    // heel strike time array : This is a list of times corresponding to heel
    // strike (0% in the gait cycle)
    // These should be set on construction and not modified.
    GainScheduledController(SimTK::Vector heel_strike_times_input,
                            std::vector<SimTK::Matrix> t_star_array_input,
                            std::vector<SimTK::Matrix> gain_array_input)
                            : OpenSim::Controller(),
                            gain_array(gain_array_input),
                            t_star_array(t_star_array_input),
                            heel_strike_times(heel_strike_times_input)
    {
    }

    void computeControls( const SimTK::State& s, SimTK::Vector& controls) const
    {
        double t = s.getTime();

        // Pointers to the joint torque actuators

        const OpenSim::CoordinateActuator* TB = dynamic_cast<const OpenSim::CoordinateActuator*> ( &getActuatorSet().get("TB") );
        const OpenSim::CoordinateActuator* TC = dynamic_cast<const OpenSim::CoordinateActuator*> ( &getActuatorSet().get("TC") );
        const OpenSim::CoordinateActuator* TD = dynamic_cast<const OpenSim::CoordinateActuator*> ( &getActuatorSet().get("TD") );
        const OpenSim::CoordinateActuator* TE = dynamic_cast<const OpenSim::CoordinateActuator*> ( &getActuatorSet().get("TE") );
        const OpenSim::CoordinateActuator* TF = dynamic_cast<const OpenSim::CoordinateActuator*> ( &getActuatorSet().get("TF") );
        const OpenSim::CoordinateActuator* TG = dynamic_cast<const OpenSim::CoordinateActuator*> ( &getActuatorSet().get("TG") );

        const OpenSim::Coordinate* qb = TB->getCoordinate();
        const OpenSim::Coordinate* qc = TC->getCoordinate();
        const OpenSim::Coordinate* qd = TD->getCoordinate();
        const OpenSim::Coordinate* qe = TE->getCoordinate();
        const OpenSim::Coordinate* qf = TF->getCoordinate();
        const OpenSim::Coordinate* qg = TG->getCoordinate();

        // Build a sensor vector.
        double qb_val = qb->getValue(s);
        double qc_val = qc->getValue(s);
        double qd_val = qd->getValue(s);
        double qe_val = qe->getValue(s);
        double qf_val = qf->getValue(s);
        double qg_val = qg->getValue(s);
        double ub_val = qb->getSpeedValue(s);
        double uc_val = qc->getSpeedValue(s);
        double ud_val = qd->getSpeedValue(s);
        double ue_val = qe->getSpeedValue(s);
        double uf_val = qf->getSpeedValue(s);
        double ug_val = qg->getSpeedValue(s);

        double s_array[12] = {qb_val, qc_val, qd_val, qe_val, qf_val, qg_val,
                              ub_val, uc_val, ud_val, ue_val, uf_val, ug_val};

        SimTK::Matrix x(12, 1, *s_array);

        double percent_gait = ComputePercentGait(t, heel_strike_times);

        // Returns a q x 1 SimTK::Matrix
        SimTK::Matrix t_star = InterpolateGainArray(percent_gait, t_star_array);

        // Returns a q x p SimTK::Matrix
        SimTK::Matrix k = InterpolateGainArray(percent_gait, gain_array);

        // Compute the control values using the control law.
        SimTK::Matrix joint_torques = t_star + k * x;

        // Apply the control values to the controller.
        SimTK::Vector TB_vec(1, joint_torques(0, 0));
        SimTK::Vector TC_vec(1, joint_torques(1, 0));
        SimTK::Vector TD_vec(1, joint_torques(2, 0));
        SimTK::Vector TE_vec(1, joint_torques(3, 0));
        SimTK::Vector TF_vec(1, joint_torques(4, 0));
        SimTK::Vector TG_vec(1, joint_torques(5, 0));

        TB->addInControls(TB_vec, controls);
        TC->addInControls(TC_vec, controls);
        TD->addInControls(TD_vec, controls);
        TE->addInControls(TE_vec, controls);
        TF->addInControls(TF_vec, controls);
        TG->addInControls(TG_vec, controls);
    }

private:

    SimTK::Vector heel_strike_times;
    std::vector<SimTK::Matrix> t_star_array;
    std::vector<SimTK::Matrix> gain_array;
};


int main()
{
    // Test InterpolateGainArray

    std::vector<SimTK::Matrix> sample_gain_array(10);

    double matrix_stepper = 0.0;
    double row_stepper;
    double col_stepper;

    for (SimTK::Matrix& gain_matrix_i : sample_gain_array) {
        int num_rows = 2;
        int num_cols = 1;
        SimTK::Matrix gain_mat(num_rows, num_cols);
        matrix_stepper += 1.0;
        for (int row_idx=0; row_idx < num_rows; row_idx++){
            row_stepper = row_idx;
            for(int col_idx=0; col_idx < num_cols; col_idx++){
                col_stepper = col_idx;
                gain_mat(row_idx, col_idx) = matrix_stepper + row_stepper + col_stepper;
            }
        }
        gain_matrix_i = gain_mat;
    }

    for (const SimTK::Matrix& gain_matrix_i : sample_gain_array) {
        std::cout << gain_matrix_i << std::endl;
    }

    SimTK::Matrix answer = InterpolateGainArray(0.87, sample_gain_array);
    std::cout << answer << std::endl;

    SimTK::Vec<10, double> percent_gait_cycle;

    for (int i=0; i<10; i++){
        percent_gait_cycle(i) = i * 10.0;
    }

    std::cout << percent_gait_cycle << std::endl;

    // Test the ComputePercentGait function.

    SimTK::Vector heel_strike_times(5);
    heel_strike_times[0] = 3.25;
    heel_strike_times[1] = 4.11;
    heel_strike_times[2] = 6.77;
    heel_strike_times[3] = 8.9;
    heel_strike_times[4] = 9.8;

    std::cout << heel_strike_times << std::endl;
    std::cout << ComputePercentGait(5.0, heel_strike_times) << std::endl;
}
