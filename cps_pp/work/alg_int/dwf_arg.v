class ActionBosonArg frm_arg = {
class ActionBilinearArg bi_arg = {
FclassType fermion = F_CLASS_DWF
Array bilinears[1] = {
class BilinearDescr bilinears[0] = {
double mass =   0.1
int max_num_iter = 500
}
}
class ActionArg action_arg = {
ForceMeasure force_measure = FORCE_MEASURE_YES
string force_label = "Fermion-Dwf"
}
}
Array fermions[1] = {
class FermionDescr fermions[0] = {
double epsilon = 1.00
int chrono = 0;
double stop_rsd_md =   1.0e-06
double stop_rsd_mc =   1.0e-10
}
}
}
