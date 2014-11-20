//This contains code fragments that are no longer in use in
//the main codebase, but should rather not get lost.


////////////////////////////////////////////////////////////////////////
////            code from template class DetSDW                     ////
////                                                                ////
////////////////////////////////////////////////////////////////////////



// template<class Matrix>
// static void setReal(arma::subview<num> subv, Matrix realPart) {
//     subv = realPart;
// }
// template<class Matrix>
// static void setReal(arma::subview<cpx> subv, Matrix realPart) {
//     //warning: unnecessarily slow
//     MatCpx temp_mat = subv;
//     temp_mat.set_real(
//     subv.set_real(realPart);
// }
// template<class Matrix>
// static void setImag(arma::subview<num> subv, Matrix imagPart) {
//     subv = imagPart;
// }
// template<class Matrix>
// static void setImag(arma::subview<cpx> subv, Matrix imagPart) {
//     subv.set_imag(imagPart);
// }

