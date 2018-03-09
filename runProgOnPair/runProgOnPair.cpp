#include <string>
#include <assert.h>
#include <iostream>

// Simple utility program that multiplexes an input "x.foo" into "x.mrc" and "x.pdb",
// and runs a given program with those parameters

int main(int argc, char *argv[])
{
	int num_args = 2;
	std::cout << "hello" << std::endl;
	std::cout << argc << std::endl;
	assert(argc == num_args + 1 && "Incorrect num arguments, 2 inputs: cmd, file.foo");

	std::string program = argv[1];
	std::string input_file_name = argv[2];
	size_t period_pos = input_file_name.find_last_of('.');
	input_file_name.resize(period_pos);

	std::string mrc_file_path = input_file_name + ".mrc";
	std::string pdb_file_path = input_file_name + ".pdb";

	std::string command = program + " \"" + mrc_file_path + "\" \"" + pdb_file_path + "\"";

	system(command.c_str());
	return 0;
}

