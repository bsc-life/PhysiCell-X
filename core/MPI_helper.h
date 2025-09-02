#include<mpi.h>
#include<vector>
#include<string>
#include<map>
using namespace std;

//packing functions to avoid code replication and enhance clarity
//pack bool
void pack_buff(const bool& input, vector<char>& buffer, int& len_buffer, int& position){
    int temp_int = (input == true) ? 1 : 0;
    len_buffer = position + sizeof(int);		
	buffer.resize(len_buffer);
	MPI_Pack(&(temp_int), 1, MPI_INT, &buffer[0], len_buffer, &position, MPI_COMM_WORLD); 
}

//pack double
void pack_buff(const double& input, vector<char>& buffer, int& len_buffer, int& position){
    len_buffer = position + sizeof(double);		
	buffer.resize(len_buffer);
	MPI_Pack(&(input), 1, MPI_DOUBLE, &buffer[0], len_buffer, &position, MPI_COMM_WORLD);
}

//pack int
void pack_buff(const int& input, vector<char>& buffer, int& len_buffer, int& position){
    len_buffer = position + sizeof(int);		
	buffer.resize(len_buffer);
	MPI_Pack(&(input), 1, MPI_INT, &buffer[0], len_buffer, &position, MPI_COMM_WORLD);
}

//pack string
void pack_buff(const string& input, vector<char>& buffer, int& len_buffer, int& position){
    
	int len_str = static_cast<int>(input.length());
	len_buffer = position + sizeof(len_str) + len_str; 
	buffer.resize(len_buffer);
	MPI_Pack(&len_str, 1, MPI_INT, &buffer[0], len_buffer, &position, MPI_COMM_WORLD); 
	MPI_Pack(&input[0], len_str, MPI_CHAR, &buffer[0], len_buffer, &position, MPI_COMM_WORLD);
}

//pack vector<int>
void pack_buff(const vector<int>& input, vector<char>& buffer, int& len_buffer, int& position){
    
	int len_vector = static_cast<int>(input.size());
	len_buffer = position + sizeof(int)  + len_vector * sizeof(int);
	buffer.resize(len_buffer);
	MPI_Pack(&(len_vector), 1, MPI_INT, &buffer[0], len_buffer, &position, MPI_COMM_WORLD);
	MPI_Pack(input.data(), len_vector, MPI_INT, &buffer[0], len_buffer, &position, MPI_COMM_WORLD);
}

//pack vector<double>
void pack_buff(const vector<double>& input, vector<char>& buffer, int& len_buffer, int& position){
	int len_vector = static_cast<int>(input.size());
	len_buffer = position + sizeof(int)  + len_vector * sizeof(double);
	buffer.resize(len_buffer);
	MPI_Pack(&(len_vector), 1, MPI_INT, &buffer[0], len_buffer, &position, MPI_COMM_WORLD);
	MPI_Pack(input.data(), len_vector, MPI_DOUBLE, &buffer[0], len_buffer, &position, MPI_COMM_WORLD);
}

//Unpacking functions
//unpack bool
void unpack_buff(bool& output, const vector<char>& buffer, int len_buffer, int& position) {
    int temp_int;
    MPI_Unpack(buffer.data(), len_buffer, &position, &temp_int, 1, MPI_INT, MPI_COMM_WORLD);
    output = (temp_int != 0);
}

void unpack_buff(int& output, const vector<char>& buffer, int len_buffer, int& position) {
    MPI_Unpack(buffer.data(), len_buffer, &position, &output, 1, MPI_INT, MPI_COMM_WORLD);
}

//unpack double
void unpack_buff(double& output, const vector<char>& buffer, int len_buffer, int& position) {
    MPI_Unpack(buffer.data(), len_buffer, &position, &output, 1, MPI_DOUBLE, MPI_COMM_WORLD);
}

//unpack string
void unpack_buff(string& output, const vector<char>& buffer, int len_buffer, int& position) {
    int len_str = 0;
    MPI_Unpack(buffer.data(), len_buffer, &position, &len_str, 1, MPI_INT, MPI_COMM_WORLD);
    if (len_str > 0) {
		output.resize(len_str);
		MPI_Unpack(buffer.data(), len_buffer, &position, &output[0], len_str, MPI_CHAR, MPI_COMM_WORLD);
	}
}

//unpack vector<int>
void unpack_buff(vector<int>& output, const vector<char>& buffer, int len_buffer, int& position) {
    int len_vector = 0;
    MPI_Unpack(buffer.data(), len_buffer, &position, &len_vector, 1, MPI_INT, MPI_COMM_WORLD);
    if (len_vector > 0) {
		output.resize(len_vector);
    	MPI_Unpack(buffer.data(), len_buffer, &position, output.data(), len_vector, MPI_INT, MPI_COMM_WORLD);
	}
}

// unpack vector<double>
void unpack_buff(std::vector<double>& output, const std::vector<char>& buffer,int len_buffer, int& position) {
    int len_vector = 0;

    // Unpack the size of the vector
    MPI_Unpack(buffer.data(), buffer.size(), &position, &len_vector, 1, MPI_INT, MPI_COMM_WORLD);

	if (len_vector > 0) {
		// Resize the output vector
		output.resize(len_vector);

		// Unpack the vector's data
		MPI_Unpack(buffer.data(), buffer.size(), &position, output.data(), len_vector, MPI_DOUBLE, MPI_COMM_WORLD);
	}
}
