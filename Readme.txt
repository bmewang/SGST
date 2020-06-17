 Short Review of Matlab-functions contained in this folder. Detailed descriptions for all functions can be evoked using the Matlab command >>help. Experiments considered here only include simple shear and biaxial stretch.
 
 
1.	SurfaceStress( Model, Parameters, Exp, Stretch, Rotation) : Main function which should be used to access cauchy stresses for any model available within this library. Said Cauchy stresses are projected such that they match a given experiment which is enumerated via a list given in function "F".
2.	F(shear, Rot, experiment)	: Is responsible to create deformation gradients based on a given experiment plus orientation of fibers.
3.	MN(N, b)					: Creates a list representing the structure tensor of "N"-th order corresponding to fiber density coefficient "b".
4.	Psi4f(f,Theta)				: Unidirectional fiber strain.
5.	Rod(theta, psi, omega)		: Returns matrix which rotates about a vector specified with "theta,psi" by angle "omega".
6.	sigmaAI(N, rho, PrincipleDirection, a,b, FF)					: Cauchy Stress Tensor based on the trapezoidal rule for AI which is performed for any deformation gradient "FF".
7.	sigmaAI_raw(theta ,phi, rho, PrincipleDirection, a, b,  FF)		: Helper function for sigmaAI.
8.	sigmaAIProjected(N, rho, PrincipleDirection, a, b, Rot, Experiment,Stretch): Projects sigmaAI tensor such that the output matches a given experimental setup.
9.  sigmaHO(Parameters, FF)		: Given deformation gradient "FF" the Cauchy stress for the Holzapfel-Ogden law is returned.
10.	sigmaHOProjected(Parameters, Rot, experiment, stretch )			: Projects sigmaHO such that the output matches a given experimental setup.
11. sigmaI1(C, Theta)			: Intermediary stress factor for 1st invariant.
12. sigmaI4(f,Theta)			: Intermediary stress factor for 4th invariant.
13. sigmaI8(f,s, Theta)			: Intermediary stress factor for 8th invariant.
14. sigmaNGST(N, StructureTensorList, PrincipleDirection, a, b, F)	: Cauchy Stress Tensor for the "N"-th order GST model for any deformation gradient "FF".
15. sigmaNGSTProjected(N, StructureTensorList, PrincipleDirection, a, b, Rot, Experiment,Stretch) : Projects sigmaNGST such that the output matches a given experimental setup.
16. sigmaNSGST(N, StructureTensorList, PrincipleDirection, a, b, F) :  Cauchy Stress Tensor for the "N"-th order SGST model for any deformation gradient "FF".
17. sigmaNSGSTProjected(N, StructureTensorList, PrincipleDirection, a, b, Rot, Experiment,Stretch): Projects sigmaNSGST such that the output matches a given experimental setup.
