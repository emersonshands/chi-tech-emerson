#ifndef CHITECH_LBS_SHELL_OPERATIONS_H
#define CHITECH_LBS_SHELL_OPERATIONS_H

#include <petscksp.h>

namespace lbs
{
int MatrixAction_Ax(Mat matrix, Vec krylov_vector, Vec Av);
int WGDSA_TGDSA_PreConditionerMult(PC pc, Vec phi_input, Vec pc_output);
}//namespace lbs

#endif //CHITECH_LBS_SHELL_OPERATIONS_H
