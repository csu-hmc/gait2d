#include <OpenSim/OpenSim.h>
#include <iostream>
#include <vector>
#include <utility>
#include <stdexcept>

// class GainScheduledController : public Controller {
// OpenSim_DECLARE_CONCRETE_OBJECT(GainScheduledController, Controller);
// 
// public:
    // GainSchedhuledController(std::Vector gainArray, std::Vector tStarArray, std::Vector heelStrikeTimeArray) : gainArray( gainArray ), tStarArray( tStarArray), heelStrikeTimeVec( heelStrikeTimeVec)
    // {
// 
    // }
// 
// void computeControls( const SimTK::State& s, SimTK::Vector controls) const
// {
    // double t = s.getTime();
// 
    // double percentGait = computePercentGait(t);
// 
    // SimTK::Matrix gainMatrix 
// 
// }
// 
// void computePercentGait(double currentTime)
// {
// 
// }
// 
// std::pair<SimTK::Matrix, SimTK::Matrix> interpolateArrays(double percentGait)
// {
// 
// }
// 
// }
//


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


int main()
{
    // Test InterpolateGainArray

    std::vector<SimTK::Matrix> sample_gain_array(10);

    double matrix_stepper = 0.0;
    double row_stepper;
    double col_stepper;

    for (SimTK::Matrix& gain_matrix_i : sample_gain_array) {
        SimTK::Matrix gain_mat(2, 4);
        matrix_stepper += 1.0;
        for (int row_idx=0; row_idx < 2; row_idx++){
            row_stepper = row_idx;
            for(int col_idx=0; col_idx < 4; col_idx++){
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
