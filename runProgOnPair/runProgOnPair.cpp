#include <string>
#include <assert.h>
#include <iostream>

// Simple utility program that multiplexes an input "x.foo" into "x.ext1" and "x.ext2",
// and runs a given program with those parameters

int main(int argc, char *argv[])
{
	int num_args = 4;
	assert(argc == num_args + 1 && "Incorrect num arguments, 4 inputs: cmd, file.ext1, ext1, ext2");

	std::string program = argv[1];
	std::string input_file_name = argv[2];
	std::string ext1 = argv[3];
	std::string ext2 = argv[4];

	size_t period_pos = input_file_name.find_last_of('.');
	input_file_name.resize(period_pos);

	std::string mrc_file_path = input_file_name + "." + ext1;
	std::string pdb_file_path = input_file_name + "." + ext2;

	std::string command = program + " \"" + mrc_file_path + "\" \"" + pdb_file_path + "\"";

	system(command.c_str());
	return 0;
}

