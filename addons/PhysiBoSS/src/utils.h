#ifndef _PhysiBoSS_utils_h_
#define _PhysiBoSS_utils_h_

#include "maboss_network.h"
#include "../../../core/MPI_helper.h"
#include "../../../core/PhysiCell_basic_signaling.h"

class MaBoSSInput
{
    static const int NODE = 0;
    static const int PARAMETER = 1;

public:
    std::string physicell_name;
    int type;
    std::string intracellular_name;
    std::string intracellular_parameter;
    std::string action;
    double threshold;
    double inact_threshold;
    double scaling;
    int smoothing;
    double smoothed_value;
    bool use_for_dead;

    MaBoSSInput() = default;

    MaBoSSInput(std::string physicell_name, std::string intracellular_name, std::string action, double threshold, double inact_threshold, int smoothing, bool use_for_dead) : physicell_name(physicell_name), intracellular_name(intracellular_name), action(action), threshold(threshold), inact_threshold(inact_threshold), smoothing(smoothing), use_for_dead(use_for_dead){
        type = NODE;
        smoothed_value = 0;
    }

    MaBoSSInput(std::string physicell_name, std::string intracellular_parameter, double scaling, int smoothing, bool use_for_dead) : physicell_name(physicell_name), intracellular_parameter(intracellular_parameter), scaling(scaling), smoothing(smoothing), use_for_dead(use_for_dead) {
        type = PARAMETER;
        smoothed_value = 0;
    }

    bool isNode() { return type == NODE; }
    bool isParameter() { return type == PARAMETER; }

    void update_value(double value) {
        smoothed_value = (smoothed_value * smoothing + value)/(smoothing + 1);
    }

    bool updateNode(bool state, double value) 
    {
        double true_value;
        if (smoothing == 0) {
            true_value = value;
        } else {
            update_value(value);
            true_value = smoothed_value;
        }

        if (state) {
            if (action == "inhibition") {
                return true_value <= inact_threshold; // When the node is active, and this is an activation, the node stays true if the value is below the inact threshold

            } else {
                return true_value >= inact_threshold; // When the node is active, the node stays true if the value is above the inact threshold
            }

        } else {
            if (action == "inhibition") {
                return true_value < threshold;

            } else {
                return true_value > threshold;
            }
        }
    }

    double updateParameter(double value) {
        if (smoothing == 0) {
            return value;
        } else {
            update_value(value);
            return smoothed_value;
        }
    }

    void pack(std::vector<char>& buffer, int& len_buffer, int& position)
    {
        // pack physicell_name
        pack_buff(physicell_name, buffer, len_buffer, position);

        // pack type
        pack_buff(type, buffer, len_buffer, position);

        // pack intracellular_name
        pack_buff(intracellular_name, buffer, len_buffer, position);

        // pack intracellular_parameter
        pack_buff(intracellular_parameter, buffer, len_buffer, position);

        // pack action
        pack_buff(action, buffer, len_buffer, position);

        // pack threshold
        pack_buff(threshold, buffer, len_buffer, position);

        // pack inact_threshold
        pack_buff(inact_threshold, buffer, len_buffer, position);

        // pack scaling
        pack_buff(scaling, buffer, len_buffer, position);

        // pack smoothing
        pack_buff(smoothing, buffer, len_buffer, position);

        // pack smoothed_value
        pack_buff(smoothed_value, buffer, len_buffer, position);

        // pack use_for_dead
        pack_buff(use_for_dead, buffer, len_buffer, position);
    }

    void unpack(const std::vector<char>& buffer, int len_buffer, int& position)
    {
        // unpack physicell_name
        unpack_buff(physicell_name, buffer, len_buffer, position);

        // unpack type
        unpack_buff(type, buffer, len_buffer, position);

        // unpack intracellular_name
        unpack_buff(intracellular_name, buffer, len_buffer, position);

        // unpack intracellular_parameter
        unpack_buff(intracellular_parameter, buffer, len_buffer, position);

        // unpack action
        unpack_buff(action, buffer, len_buffer, position);

        // unpack threshold
        unpack_buff(threshold, buffer, len_buffer, position);

        // unpack inact_threshold
        unpack_buff(inact_threshold, buffer, len_buffer, position);

        // unpack scaling
        unpack_buff(scaling, buffer, len_buffer, position);

        // unpack smoothing
        unpack_buff(smoothing, buffer, len_buffer, position);

        // unpack smoothed_value
        unpack_buff(smoothed_value, buffer, len_buffer, position);

        // unpack use_for_dead
        unpack_buff(use_for_dead, buffer, len_buffer, position);
    }


    /* OG versions from Gaurov substituted above for increased modularity
    void MaBoSSInput::pack(std::vector<char>& buffer, int& len_buffer, int& position)
    {
        //pack physicell_name
        int len_str = physicell_name.length();
        len_buffer = position + sizeof(len_str) + len_str;
        buffer.resize(len_buffer);
        MPI_Pack(&len_str, 1, MPI_INT, &buffer[0], len_buffer, &position, MPI_COMM_WORLD);
        MPI_Pack(&physicell_name[0], len_str, MPI_CHAR, &buffer[0], len_buffer, &position, MPI_COMM_WORLD);

        //pack type 
        len_buffer = position + sizeof(int);
        buffer.resize(len_buffer);
        int temp_int = type;
        MPI_Pack(&temp_int, 1, MPI_INT, &buffer[0], len_buffer, &position, MPI_COMM_WORLD);

        //pack intracellular_name
        int len_str = intracellular_name.length();
        len_buffer = position + sizeof(len_str) + len_str;
        buffer.resize(len_buffer);
        MPI_Pack(&len_str, 1, MPI_INT, &buffer[0], len_buffer, &position, MPI_COMM_WORLD);
        MPI_Pack(&intracellular_name[0], len_str, MPI_CHAR, &buffer[0], len_buffer, &position, MPI_COMM_WORLD);

        //pack intracellular_parameter
        int len_str = intracellular_parameter.length();
        len_buffer = position + sizeof(len_str) + len_str;
        buffer.resize(len_buffer);
        MPI_Pack(&len_str, 1, MPI_INT, &buffer[0], len_buffer, &position, MPI_COMM_WORLD);
        MPI_Pack(&intracellular_parameter[0], len_str, MPI_CHAR, &buffer[0], len_buffer, &position, MPI_COMM_WORLD);

        //pack action
        //pack physicell_name
        int len_str = action.length();
        len_buffer = position + sizeof(len_str) + len_str;
        buffer.resize(len_buffer);
        MPI_Pack(&len_str, 1, MPI_INT, &buffer[0], len_buffer, &position, MPI_COMM_WORLD);
        MPI_Pack(&action[0], len_str, MPI_CHAR, &buffer[0], len_buffer, &position, MPI_COMM_WORLD);

        //pack threshold
        len_buffer = position + sizeof(double);
        buffer.resize(len_buffer);
        MPI_Pack(&threshold, 1, MPI_DOUBLE, &buffer[0], len_buffer, &position, MPI_COMM_WORLD);

        //pack inact_threshold
        len_buffer = position + sizeof(double);
        buffer.resize(len_buffer);
        MPI_Pack(&inact_threshold, 1, MPI_DOUBLE, &buffer[0], len_buffer, &position, MPI_COMM_WORLD);

        //pack scaling
        len_buffer = position + sizeof(double);
        buffer.resize(len_buffer);
        MPI_Pack(&scaling, 1, MPI_DOUBLE, &buffer[0], len_buffer, &position, MPI_COMM_WORLD);

        //pack smoothing
        len_buffer = position + sizeof(int);
        buffer.resize(len_buffer);
        MPI_Pack(&smoothing, 1, MPI_INT, &buffer[0], len_buffer, &position, MPI_COMM_WORLD);

        //pack smoothed_value
        len_buffer = position + sizeof(double);
        buffer.resize(len_buffer);
        MPI_Pack(&smoothed_value, 1, MPI_DOUBLE, &buffer[0], len_buffer, &position, MPI_COMM_WORLD);

        //pack use_for_dead
        int use_for_dead_int = (use_for_dead ? 1 : 0);
        len_buffer = position + sizeof(int);
        buffer.resize(len_buffer);
        MPI_Pack(&use_for_dead_int, 1, MPI_INT, &buffer[0], len_buffer, &position, MPI_COMM_WORLD);

    }

    void MaBoSSInput::unpack(const std::vector<char>& buffer, int& position)
    {
        // unpack physicell_name
        int len_str;
        MPI_Unpack(buffer.data(), buffer.size(), &position, &len_str, 1, MPI_INT, MPI_COMM_WORLD);
        physicell_name.resize(len_str);
        MPI_Unpack(buffer.data(), buffer.size(), &position, &physicell_name[0], len_str, MPI_CHAR, MPI_COMM_WORLD);

        // unpack type
        int temp_int;
        MPI_Unpack(buffer.data(), buffer.size(), &position, &temp_int, 1, MPI_INT, MPI_COMM_WORLD);
        type = temp_int;

        // unpack intracellular_name
        MPI_Unpack(buffer.data(), buffer.size(), &position, &len_str, 1, MPI_INT, MPI_COMM_WORLD);
        intracellular_name.resize(len_str);
        MPI_Unpack(buffer.data(), buffer.size(), &position, &intracellular_name[0], len_str, MPI_CHAR, MPI_COMM_WORLD);

        // unpack intracellular_parameter
        MPI_Unpack(buffer.data(), buffer.size(), &position, &len_str, 1, MPI_INT, MPI_COMM_WORLD);
        intracellular_parameter.resize(len_str);
        MPI_Unpack(buffer.data(), buffer.size(), &position, &intracellular_parameter[0], len_str, MPI_CHAR, MPI_COMM_WORLD);

        // unpack action
        MPI_Unpack(buffer.data(), buffer.size(), &position, &len_str, 1, MPI_INT, MPI_COMM_WORLD);
        action.resize(len_str);
        MPI_Unpack(buffer.data(), buffer.size(), &position, &action[0], len_str, MPI_CHAR, MPI_COMM_WORLD);

        // unpack threshold
        MPI_Unpack(buffer.data(), buffer.size(), &position, &threshold, 1, MPI_DOUBLE, MPI_COMM_WORLD);

        // unpack inact_threshold
        MPI_Unpack(buffer.data(), buffer.size(), &position, &inact_threshold, 1, MPI_DOUBLE, MPI_COMM_WORLD);

        // unpack scaling
        MPI_Unpack(buffer.data(), buffer.size(), &position, &scaling, 1, MPI_DOUBLE, MPI_COMM_WORLD);

        // unpack smoothing
        MPI_Unpack(buffer.data(), buffer.size(), &position, &smoothing, 1, MPI_INT, MPI_COMM_WORLD);

        // unpack smoothed_value
        MPI_Unpack(buffer.data(), buffer.size(), &position, &smoothed_value, 1, MPI_DOUBLE, MPI_COMM_WORLD);

        // unpack use_for_dead
        int use_for_dead_int;
        MPI_Unpack(buffer.data(), buffer.size(), &position, &use_for_dead_int, 1, MPI_INT, MPI_COMM_WORLD);
        use_for_dead = (use_for_dead_int != 0);
    }*/
};

class MaBoSSOutput
{
public:
    std::string physicell_name;
    std::string intracellular_name;
    std::string action;
    double value;
    double base_value;
    int smoothing;
    double probability;
    bool initialized = false;
    int steepness;
    bool use_for_dead;

    MaBoSSOutput() = default;

    MaBoSSOutput(std::string physicell_name, std::string intracellular_name,
                std::string action, double value, double base_value,
                int smoothing, int steepness, bool use_for_dead)
        : physicell_name(physicell_name), intracellular_name(intracellular_name),
        action(action), value(value), base_value(base_value),
        smoothing(smoothing), steepness(steepness), use_for_dead(use_for_dead) {
    probability = 0.5;
    }

    double update_probability(bool test) {
    if (!initialized) {
        probability = test;
        initialized = true;
    } else
        probability =
            (probability * smoothing + (test ? 1.0 : 0.0)) / (smoothing + 1.0);

    return probability;
    }

    double update(bool test) {

    double hill_input;

    if (smoothing == 0) {
        hill_input = test ? 1.0 : 0.0;
    } else {
        hill_input = update_probability(test);
    }

    if (action == "activation") {
        double hill = PhysiCell::Hill_response_function(hill_input * 2, 1, steepness);
        return (value - base_value) * hill + base_value;
    } else if (action == "inhibition") {
        double hill = PhysiCell::Hill_response_function(hill_input * 2, 1, steepness);
        return ((value - base_value) * (1 - hill)) + base_value;
    }

    return base_value;
    }

    void pack(std::vector<char>& buffer, int& len_buffer, int& position)
    {
        // pack physicell_name
        pack_buff(physicell_name, buffer, len_buffer, position);

        // pack intracellular_name
        pack_buff(intracellular_name, buffer, len_buffer, position);

        // pack action
        pack_buff(action, buffer, len_buffer, position);

        // pack value
        pack_buff(value, buffer, len_buffer, position);

        // pack base_value
        pack_buff(base_value, buffer, len_buffer, position);

        // pack smoothing
        pack_buff(smoothing, buffer, len_buffer, position);

        // pack probability
        pack_buff(probability, buffer, len_buffer, position);

        // pack initialized
        pack_buff(initialized, buffer, len_buffer, position);

        // pack steepness
        pack_buff(steepness, buffer, len_buffer, position);

        // pack use_for_dead
        pack_buff(use_for_dead, buffer, len_buffer, position);
    }

    void unpack(const std::vector<char>& buffer, int len_buffer, int& position)
    {
        // unpack physicell_name
        unpack_buff(physicell_name, buffer, len_buffer, position);

        // unpack intracellular_name
        unpack_buff(intracellular_name, buffer, len_buffer, position);

        // unpack action
        unpack_buff(action, buffer, len_buffer, position);

        // unpack value
        unpack_buff(value, buffer, len_buffer, position);

        // unpack base_value
        unpack_buff(base_value, buffer, len_buffer, position);

        // unpack smoothing
        unpack_buff(smoothing, buffer, len_buffer, position);

        // unpack probability
        unpack_buff(probability, buffer, len_buffer, position);

        // unpack initialized
        unpack_buff(initialized, buffer, len_buffer, position);

        // unpack steepness
        unpack_buff(steepness, buffer, len_buffer, position);

        // unpack use_for_dead
        unpack_buff(use_for_dead, buffer, len_buffer, position);
    }

    /*
    void MaBoSSOutput::pack(std::vector<char>& buffer, int& len_buffer, int& position)
    {
        //pack physicell_name
        int len_str = physicell_name.length();
        len_buffer = position + sizeof(len_str) + len_str;
        buffer.resize(len_buffer);
        MPI_Pack(&len_str, 1, MPI_INT, &buffer[0], len_buffer, &position, MPI_COMM_WORLD);
        MPI_Pack(&physicell_name[0], len_str, MPI_CHAR, &buffer[0], len_buffer, &position, MPI_COMM_WORLD);
        
        //pack intracellular_name
        len_str = intracellular_name.length();
        len_buffer = position + sizeof(len_str) + len_str;
        buffer.resize(len_buffer);
        MPI_Pack(&len_str, 1, MPI_INT, &buffer[0], len_buffer, &position, MPI_COMM_WORLD);
        MPI_Pack(&intracellular_name[0], len_str, MPI_CHAR, &buffer[0], len_buffer, &position, MPI_COMM_WORLD);

        //pack action
        len_str = action.length();
        len_buffer = position + sizeof(len_str) + len_str;
        buffer.resize(len_buffer);
        MPI_Pack(&len_str, 1, MPI_INT, &buffer[0], len_buffer, &position, MPI_COMM_WORLD);
        MPI_Pack(&action[0], len_str, MPI_CHAR, &buffer[0], len_buffer, &position, MPI_COMM_WORLD);

        //pack value 
        len_buffer = position + sizeof(double);
        buffer.resize(len_buffer);
        MPI_Pack(&value, 1, MPI_DOUBLE, &buffer[0], len_buffer, &position, MPI_COMM_WORLD);

        //pack base value 
        len_buffer = position + sizeof(double);
        buffer.resize(len_buffer);
        MPI_Pack(&base_value, 1, MPI_DOUBLE, &buffer[0], len_buffer, &position, MPI_COMM_WORLD);

        //pack smoothing
        len_buffer = position + sizeof(int);
        buffer.resize(len_buffer);
        MPI_Pack(&smoothing, 1, MPI_INT, &buffer[0], len_buffer, &position, MPI_COMM_WORLD);

        //pack probability
        len_buffer = position + sizeof(double);
        buffer.resize(len_buffer);
        MPI_Pack(&probability, 1, MPI_DOUBLE, &buffer[0], len_buffer, &position, MPI_COMM_WORLD);

        //pack initialized
        int initialized_int = (initialized ? 1 : 0);
        len_buffer = position + sizeof(int);
        buffer.resize(len_buffer);
        MPI_Pack(&initialized_int, 1, MPI_INT, &buffer[0], len_buffer, &position, MPI_COMM_WORLD);

        //pack steepness
        len_buffer = position + sizeof(int);
        buffer.resize(len_buffer);
        MPI_Pack(&steepness, 1, MPI_INT, &buffer[0], len_buffer, &position, MPI_COMM_WORLD);

        //pack use_for_dead
        int use_for_dead_int = (use_for_dead ? 1 : 0);
        len_buffer = position + sizeof(int);
        buffer.resize(len_buffer);
        MPI_Pack(&use_for_dead_int, 1, MPI_INT, &buffer[0], len_buffer, &position, MPI_COMM_WORLD);
    }

    void MaBoSSOutput::unpack(const std::vector<char>& buffer, int& position)
    {
        int len_str;

        // Unpack physicell_name
        MPI_Unpack(buffer.data(), buffer.size(), &position, &len_str, 1, MPI_INT, MPI_COMM_WORLD);
        physicell_name.resize(len_str);
        MPI_Unpack(buffer.data(), buffer.size(), &position, &physicell_name[0], len_str, MPI_CHAR, MPI_COMM_WORLD);

        // Unpack intracellular_name
        MPI_Unpack(buffer.data(), buffer.size(), &position, &len_str, 1, MPI_INT, MPI_COMM_WORLD);
        intracellular_name.resize(len_str);
        MPI_Unpack(buffer.data(), buffer.size(), &position, &intracellular_name[0], len_str, MPI_CHAR, MPI_COMM_WORLD);

        // Unpack action
        MPI_Unpack(buffer.data(), buffer.size(), &position, &len_str, 1, MPI_INT, MPI_COMM_WORLD);
        action.resize(len_str);
        MPI_Unpack(buffer.data(), buffer.size(), &position, &action[0], len_str, MPI_CHAR, MPI_COMM_WORLD);

        // Unpack value
        MPI_Unpack(buffer.data(), buffer.size(), &position, &value, 1, MPI_DOUBLE, MPI_COMM_WORLD);

        // Unpack base_value
        MPI_Unpack(buffer.data(), buffer.size(), &position, &base_value, 1, MPI_DOUBLE, MPI_COMM_WORLD);

        // Unpack smoothing
        MPI_Unpack(buffer.data(), buffer.size(), &position, &smoothing, 1, MPI_INT, MPI_COMM_WORLD);

        // Unpack probability
        MPI_Unpack(buffer.data(), buffer.size(), &position, &probability, 1, MPI_DOUBLE, MPI_COMM_WORLD);

        // Unpack initialized
        int initialized_int;
        MPI_Unpack(buffer.data(), buffer.size(), &position, &initialized_int, 1, MPI_INT, MPI_COMM_WORLD);
        initialized = (initialized_int != 0);

        // Unpack steepness
        MPI_Unpack(buffer.data(), buffer.size(), &position, &steepness, 1, MPI_INT, MPI_COMM_WORLD);

        // Unpack use_for_dead
        int use_for_dead_int;
        MPI_Unpack(buffer.data(), buffer.size(), &position, &use_for_dead_int, 1, MPI_INT, MPI_COMM_WORLD);
        use_for_dead = (use_for_dead_int != 0);
    } */

};
#endif
