//------------------------------------------------------------------------
// Code ported from:
// Johan Clausen
// Esbjerg Department of Engineering
// Aalborg University
// July 2005
//------------------------------------------------------------------------


using ISAAR.MSolve.Materials.Interfaces;
using ISAAR.MSolve.Numerical.LinearAlgebra;
using ISAAR.MSolve.Numerical.LinearAlgebra.Interfaces;
using System;


namespace ISAAR.MSolve.FEM.Materials
{
    public class MohrCoulombMaterial : IIsotropicContinuumMaterial3D
    {
        private const double PoissonRatioForIncompressibleSolid = 0.5;
        private bool modified;
        //private readonly double[,] elasticConstitutiveMatrix;
        private int plasticRegion;
        private double[,] constitutiveMatrix = new double[6, 6];
        private double[,] constitutiveMatrixNew = new double[6, 6];
        private double[] incrementalStrains = new double[6];
        private double[] stresses = new double[6];
        private double[] stressesNew = new double[6];
        private readonly double youngModulus, shearModulus, poissonRatio, cohesion, friction, dilation;

        public StressStrainVectorContinuum3D Stresses
        {
            get { return new StressStrainVectorContinuum3D(stressesNew); }
        }

        public ElasticityTensorContinuum3D ConstitutiveMatrix
        {
            get { return new ElasticityTensorContinuum3D(constitutiveMatrix); }
        }

        public void UpdateMaterial(StressStrainVectorContinuum3D strainsIncrement)
        {
            Array.Copy(strainsIncrement.Data, this.incrementalStrains, 6);
            this.CalculateNextStressStrainPoint();
        }

        public void ClearState()
        {
            constitutiveMatrix = new double[6, 6];
            constitutiveMatrixNew = new double[6, 6];
            incrementalStrains = new double[6];
            stresses = new double[6];
            stressesNew = new double[6];

            var Dinv = new double[6, 6];
            DlinElas(youngModulus, poissonRatio, 6, constitutiveMatrix, Dinv);
        }

        public void SaveState()
        {
            Array.Copy(this.constitutiveMatrixNew, this.constitutiveMatrix, 6 * 6);
            Array.Copy(this.stressesNew, this.stresses, 6);
        }

        public void ClearStresses()
        {
            Array.Clear(this.stresses, 0, 6);
            Array.Clear(this.stressesNew, 0, 6);
        }

        public int ID
        {
            get { return 998; }
        }

        public bool Modified
        {
            get { return modified; }
        }

        public void ResetModified()
        {
            modified = false;
        }

        public double YoungModulus
        {
            get { return youngModulus; }
            set { throw new InvalidOperationException(); }
        }

        public double PoissonRatio
        {
            get { return poissonRatio; }
            set { throw new InvalidOperationException(); }
        }

        public double[] Coordinates { get; set; }
        public double Cohesion { get { return cohesion; } }
        public double Friction { get { return friction; } }
        public double Dilation { get { return dilation; } }

        public object Clone()
        {
            var constitutiveMatrixCopy = new double[6, 6];
            Array.Copy(constitutiveMatrix, constitutiveMatrixCopy, 36);
            var strainsCopy = new double[incrementalStrains.Length];
            Array.Copy(incrementalStrains, strainsCopy, incrementalStrains.Length);
            var stressesCopy = new double[stresses.Length];
            Array.Copy(stresses, stressesCopy, stresses.Length);

            var m = new MohrCoulombMaterial(this.youngModulus, this.poissonRatio, this.cohesion, this.friction, this.dilation)
            {
                modified = this.Modified,
                constitutiveMatrix = constitutiveMatrixCopy,
                incrementalStrains = strainsCopy,
                stresses = stressesCopy
            };
            return m;
        }

        public MohrCoulombMaterial(double youngModulus, double poissonRatio, double cohesion, double friction, double dilation)
        {
            this.youngModulus = youngModulus;

            if (poissonRatio == PoissonRatioForIncompressibleSolid)
            {
                throw new ArgumentException(
                    "Poisson ratio cannot be" + PoissonRatioForIncompressibleSolid + "(incompressible solid)");
            }

            this.poissonRatio = poissonRatio;
            this.cohesion = cohesion;
            this.friction = friction;
            this.dilation = dilation;

            this.shearModulus = this.YoungModulus / (2 * (1 + this.PoissonRatio));
            var Dinv = new double[6, 6];
            DlinElas(youngModulus, poissonRatio, 6, constitutiveMatrix, Dinv);

        }

        private void CalculateNextStressStrainPoint()
        {
            var E = youngModulus;
            var nu = poissonRatio;
            var coh = cohesion;
            var phi_d = friction; // degrees
            var psi_d = dilation; // degrees
            var phi = Math.PI * phi_d / 180d;
            var psi = Math.PI * psi_d / 180d;
            const int NTENS = 6;
            var D = new double[6, 6];
            var Dinv = new double[6, 6];
            DlinElas(E, nu, NTENS, D, Dinv);

            var SigB = new double[6];
            for (int i = 0; i < 6; i++)
            {
                SigB[i] = this.stresses[i];
                for (int j = 0; j < 6; j++)
                    SigB[i] += D[i, j] * this.incrementalStrains[j];
            }

            //-----------------------------------------------------------------------
            //     Stress update and creation of consistent constitutive matrix
            //-----------------------------------------------------------------------

            var PlastPar = new double[3];
            PlastPar[0] = (1d + Math.Sin(phi)) / (1d - Math.Sin(phi));
            PlastPar[1] = 2d * coh * Math.Sqrt(PlastPar[0]);
            PlastPar[2] = (1d + Math.Sin(psi)) / (1d - Math.Sin(psi));

            MohrCoulombStressReturn(SigB, NTENS, PlastPar, D, Dinv, stressesNew, constitutiveMatrixNew, out plasticRegion);
            this.modified = CheckIfConstitutiveMatrixChanged();
        }

        private bool CheckIfConstitutiveMatrixChanged()
        {
            for (int i = 0; i < 6; i++)
                for (int j = 0; j < 6; j++)
                    if (Math.Abs(constitutiveMatrix[i, j] - constitutiveMatrixNew[i, j]) > 1e-10)
                        return true;

            return false;
        }

        private void DlinElas(double E, double nu, int nsigma, double[,] D, double[,] Dinv)
        {
            //----- Plane stress with three stress components -----------------------
            if (nsigma == 3)
            {
                var coeff = E / (1 - nu * nu);
                D[1, 1] = coeff;
                D[1, 2] = coeff * nu;
                D[2, 1] = coeff * nu;
                D[2, 2] = coeff;
                D[3, 3] = coeff * (1 - nu) * 0.5;
                Dinv[1, 1] = 1d / E;
                Dinv[1, 2] = -nu / E;
                Dinv[2, 1] = -nu / E;
                Dinv[2, 2] = 1d / E;
                Dinv[3, 3] = 2d * (1 + nu) / E;
            }
            //----- Plane stress, strain or axisymmetry (four stress components) ---
            if (nsigma == 4)
            {
                var coeff = E / ((1 + nu) * (1 - 2d * nu));
                D[0, 0] = coeff * (1d - nu);
                D[0, 1] = coeff * nu;
                D[0, 2] = coeff * nu;
                D[1, 0] = coeff * nu;
                D[1, 1] = coeff * (1d - nu);
                D[1, 2] = coeff * nu;
                D[2, 0] = coeff * nu;
                D[2, 1] = coeff * nu;
                D[2, 2] = coeff * (1d - nu);
                D[3, 3] = coeff * (1d - 2d * nu) * 0.5;

                Dinv[0, 0] = 1d / E;
                Dinv[0, 1] = -nu / E;
                Dinv[0, 2] = -nu / E;
                Dinv[1, 0] = -nu / E;
                Dinv[1, 1] = 1d / E;
                Dinv[1, 2] = -nu / E;
                Dinv[2, 0] = -nu / E;
                Dinv[2, 1] = -nu / E;
                Dinv[2, 2] = 1d / E;
                Dinv[3, 3] = 2d * (1 + nu) / E;
            }
            if (nsigma == 6)
            {
                var c = E / ((1 + nu) * (1 - 2d * nu));
                var g = E / (2d * (1 + nu));
                D[0, 0] = (1 - nu) * c;
                D[0, 1] = nu * c;
                D[0, 2] = D[0, 1];
                D[1, 0] = D[0, 1];
                D[1, 1] = D[0, 0];
                D[1, 2] = D[0, 1];
                D[2, 0] = D[0, 1];
                D[2, 1] = D[0, 1];
                D[2, 2] = D[0, 0];
                D[3, 3] = g;
                D[4, 4] = g;
                D[5, 5] = g;

                Dinv[0, 0] = 1d / E;
                Dinv[0, 1] = -nu / E;
                Dinv[0, 2] = -nu / E;
                Dinv[1, 0] = -nu / E;
                Dinv[1, 1] = 1d / E;
                Dinv[1, 2] = -nu / E;
                Dinv[2, 0] = -nu / E;
                Dinv[2, 1] = -nu / E;
                Dinv[2, 2] = 1d / E;
                Dinv[3, 3] = 2d * (1 + nu) / E;
                Dinv[4, 4] = 2d * (1 + nu) / E;
                Dinv[5, 5] = 2d * (1 + nu) / E;
            }
        }

        //------------------------------------------------------------------------
        // Wrapping subroutine for return mapping with a linearly elastic perfectly
        // plastic Mohr-Coulomb material. Performed in principal coordinates.
        //
        // The Mohr-Coulomb criterion is defined by
        //		f = k*sigP(1) - sigP(3) - comp = 0
        //
        // and the corresponding plastic potential
        //       g = m*sigP(1) - sigP(3)
        // 
        //
        // INPUT
        //  - Name -		-type and description --
        //	Sigma		real array(nsigma). Stress vector. Size is dependent on the type of stress.
        //				  The ordering of the components must be as follows:
        //				    Plane situations (plane stress, plane strain and axisymmetry):
        //					Sigma = [sig_x sig_y sig_z tau_xy] (i.e. sig_z is out-of-plane)
        //					General 3D:  Sigma = [sig_x sig_y sig_z tau_xy tau_xz tau_yz]
        //	nsigma		integer,scalar. Number of stress components, i.e. the size of Sigma
        //	PlasPar		real, array(3). Vector of plastic material parameters [k comp m]
        //							k: Describes internal friction, k = (1+sin(phi))/(1-sin(phi))
        //							comp: Uniaxial compressive strength. comp = 2*c*sqrt(k)
        //							m: Describes dilation. m = (1+sin(psi))/(1-sin(psi))
        //	D			real, array(nsigma,nsigma) Elastic constitutive matrix.
        //	Dinv		real, array(nsigma,nsigma) Inverted elastic constitutive matrix (compliance matrix).
        //
        // OUTPUT
        //	Sigma_up	real, array(nsigma) Stress vector containg the updated stresses.
        //	Depc		real, array(nsigma,nsigma) Consistent constitutive elasto-plastic matrix
        //				  If the material is not yielding Depc = D
        //	region		integer, scalar. Number of the region of the returned stress
        //					region = 0: the elastic region, i.e. no yielding
        //					region = 1: Return to the Mohr-Coulomb surface
        //					region = 2: Return to the triaxial compressive line, sigp_1 = sigp_2
        //					region = 3: Return to the triaxial tensile line, sigp_2 = sigp_3
        //					region = 4: Return to the apex, sigp_1 = sigp_2 = sigp_3.
        //
        //------------------------------------------------------------------------
        // Johan Clausen
        // Department of Civil Engineering
        // Aalborg University
        // May 2008
        //------------------------------------------------------------------------
        private void MohrCoulombStressReturn(double[] Sigma, int nsigma, double[] PlasPar, double[,] D, double[,] Dinv, double[] Sigma_up,
            double[,] Depc, out int region)
        {
            int ouplP, s1, s2, nshear;
            double f;
            double[,] psi = new double[3, 3];
            double[,] DepP = new double[nsigma, nsigma];
            double[,] DepcP = new double[nsigma, nsigma];
            double[,] A = new double[nsigma, nsigma];
            double[,] Atrans = new double[nsigma, nsigma];
            double[,] MatProd = new double[nsigma, nsigma];
            double[,] T = new double[nsigma, nsigma];
            double[,] Tshear = new double[nsigma - 3, nsigma - 3];
            double[,] Tshearinv = new double[nsigma - 3, nsigma - 3];
            double[] SigP = new double[3];
            double[] SigP_up = new double[3];
            double[] Fnorm;
            double[] Gnorm;
            double[] Lfdir;
            double[] Lgdir;
            double[] SiPla = new double[3];

            //--- Settings regarding calculation of the constitutive matrix, DepP -------------
            const int DepFormLine = 1; // Form of the consistent constitutive matrix on a line
                                       //! DepFormLine = 0: Double singular constitutive matrix
                                       //! DepFormLine = 1: Modified "almost double singular" matrix on a line ("Tiln�_b=axr" in the MatLab code)
                                       //! DepFormLine = 2: A single singular matrix in the Koiter direction, i.e. the direction of the plastic strain.
            const double beta = 20d; // The value used in the modification of DepP when DepFormLine = 1
            const int DepFormApex = 1; // Form of the consistent constitutive matrix on the apex
                                       //! DepFormApex = 0: DepP is the zero matrix
                                       //! DepFormApex = 1: A modified single singular DepP in the Koiter direction ("DepKo/fak" in the MatLab code)
                                       //! DepFormApex = 2: A modified double singular DepP in the Koiter and the hydrostatic direction ("DepKoiteP" in the MatLab code)
            const double alpha = 100d; // Factor in the modified apex-Depc-formulation. Used when DepFormApex = 1 or DepFormApex = 2

            PrinStressAna(Sigma, nsigma, SigP); // Principal stresses SigP
                                                //	--- Value of Yield function f ---
                                                //!	PlasPar = [k comp m]
            f = PlasPar[0] * SigP[0] - SigP[2] - PlasPar[1]; // f = k*SigP_1 - SigP_3 - comp
            s1 = 0;
            s2 = 0;
            ouplP = 0;
            if (f > 0d)
            {
                //		--- Position of out-of plane or hoop stress in principal stress vector SigP ----
                if (nsigma == 4)
                {
                    if (Sigma[2] == SigP[0])
                    {
                        ouplP = 1; s1 = 2; s2 = 3;
                    }
                    else if (Sigma[2] == SigP[1])
                    {
                        ouplP = 2; s1 = 1; s2 = 3;
                    }
                    else if (Sigma[2] == SigP[2])
                    {
                        ouplP = 3; s1 = 1; s2 = 2;
                    }
                }

                //		--- Stress return in principal coordinates --
                PrinRetMoCo(SigP, f, PlasPar, D, nsigma, SigP_up, out region);
                for (int i = 0; i < 3; i++)
                    SiPla[i] = SigP[i] - SigP_up[i]; // Plastic corrector stress in principal stress space

                //		--- Modification matrix, T ----------------------
                nshear = nsigma - 3; // Number of shear components
                TshearPrinPerfect(SigP, SigP_up, nshear, s1 - 1, s2 - 1, Tshear, Tshearinv); // Forms the shear part modification matrix T
                T[0, 0] = 1;
                T[1, 1] = 1;
                T[2, 2] = 1;
                for (int i = 3; i < nsigma; i++)
                    for (int j = 3; j < nsigma; j++)
                        T[i, j] = Tshear[i - 3, j - 3];

                //		--- Infinitesimal constitutive matrix, DepP -----------------	
                if (region == 1) // Return to the yield surface 
                {
                    for (int i = 3; i < nsigma; i++)
                        for (int j = 3; j < nsigma; j++)
                            DepP[i, j] = D[i, j];     // Shear components of consistent constitutive matrix in principal stress space

                    // PlasPar = [k comp m]
                    Fnorm = new double[] { PlasPar[0], 0d, -1d }; // Yield plane normal
                    Gnorm = new double[] { PlasPar[2], 0d, -1d }; // Potential plane normal
                    double[,] DCopy = new double[3, 3];
                    double[,] DepPCopy = new double[3, 3];
                    for (int i = 0; i < 3; i++)
                        for (int j = 0; j < 3; j++)
                        {
                            DCopy[i, j] = D[i, j];
                            DepPCopy[i, j] = DepP[i, j];
                        }
                    FormDepPerfect(DCopy, Fnorm, Gnorm, 3, DepPCopy); // Normal components of consistent constitutive matrix in principal stress space
                    for (int i = 0; i < 3; i++)
                        for (int j = 0; j < 3; j++)
                            DepP[i, j] = DepPCopy[i, j];
                }
                else if (region == 2) // return to a line 1, triaxial compression, sigp1 = sigp2
                {
                    Lfdir = new double[] { 1d, 1d, PlasPar[0] }; // Edge line direction
                    Lgdir = new double[] { 1d, 1d, PlasPar[2] }; // Potential edge line direction
                    double[,] DinvCopy = new double[3, 3];
                    for (int i = 0; i < 3; i++)
                        for (int j = 0; j < 3; j++)
                            DinvCopy[i, j] = Dinv[i, j];
                    FormModDepLine(DepFormLine, beta, SiPla, Lfdir, Lgdir, DinvCopy, D, Dinv, nsigma, DepP);
                }
                else if (region == 3) // return to a line 1, triaxial extension, sigp2 = sigp3
                {
                    Lfdir = new double[] { 1d, PlasPar[0], PlasPar[0] }; // Edge line direction
                    Lgdir = new double[] { 1d, PlasPar[2], PlasPar[2] }; // Potential edge line direction
                    double[,] DinvCopy = new double[3, 3];
                    for (int i = 0; i < 3; i++)
                        for (int j = 0; j < 3; j++)
                            DinvCopy[i, j] = Dinv[i, j];
                    FormModDepLine(DepFormLine, beta, SiPla, Lfdir, Lgdir, DinvCopy, D, Dinv, nsigma, DepP);
                }
                else // apex return
                {
                    double[,] DinvCopy = new double[3, 3];
                    for (int i = 0; i < 3; i++)
                        for (int j = 0; j < 3; j++)
                            DinvCopy[i, j] = Dinv[i, j];
                    FormModDepApex(DepFormApex, alpha, SiPla, DinvCopy, D, Dinv, nsigma, DepP);
                }

                //		--- Consistent constitutive matrix, DepcP -------	
                matmul(T, DepP, DepcP);
                //		--- Tranformation matrix A ------------------
                PrinDirect(Sigma, nsigma, SigP, psi);
                TransMatrix(psi, nsigma, ouplP, A);
                //		--- Coordinate transformation ---------------
                for (int i = 0; i < nsigma; i++)
                    for (int j = 0; j < nsigma; j++)
                        Atrans[i, j] = A[j, i];

                Array.Clear(Sigma_up, 0, nsigma);
                double[,] AtransCopy = new double[nsigma, 3];
                for (int i = 0; i < nsigma; i++)
                    for (int j = 0; j < 3; j++)
                        AtransCopy[i, j] = Atrans[i, j];
                matmul(AtransCopy, SigP_up, Sigma_up);
                matmul(Atrans, DepcP, MatProd);
                matmul(MatProd, A, Depc);
            }
            else // no yielding
            {
                for (int i = 0; i < nsigma; i++)
                {
                    Sigma_up[i] = Sigma[i];
                    for (int j = 0; j < nsigma; j++)
                        Depc[i, j] = D[i, j];
                }
                region = 0;
            }
        }


        //-----------------------------------------------------------------------
        // Stress return for Mohr-Coulomb plasticity in principal stress space.
        // Includes non-associated plasticity, but no hardening is allowed.
        // The Mohr-Coulomb criterion is written as f = k*sigP_1 - sigP_3 - comp = 0
        // The plastic potential is m = m*sigP_1 - sigP_3
        //
        // INPUT
        //  - Name - -type,size -  -- Description --
        //     SigP  (real,ar(3))  Principal predictor stresses in descending order
        //                   SigP = [sigP_1 sigP_2 sigP_3]
        //     f	  (real,sc)   Value of the yield function given by
        //                   f = k*sigP_1 - sigP_3 - comp
        //     PlasPar (real,ar(3))  Vector of plastic material parameters  = [k comp m]
        //                   k: Describes internal friction, k = (1+sin(phi))/(1-sin(phi))
        //                   comp: Uniaxial compressive strength. comp = 2*c*sqrt(k)
        //                   m: Describes dilation. m = (1+sin(psi))/(1-sin(psi))
        //                   psi = Friction angle, c = cohesion, psi = dilation angle
        //
        //     Dfull (real,ar(nsigma,nsigma)) Elastic constitutive matrix
        //     nsigma  (integer,sc)  Number of stress components, Dfull
        //
        // OUTPUT
        //     SigP_up (real,ar(3))  Updated principal stresses, i.e. they obey the
        //                   Mohr-Coulomb yield criterion
        //     region  (integer,sc)  Type of return, i.e. the region of the
        //                   predictor stress
        //
        //----------------------------------------------------------------------
        // Johan Clausen
        // Esbjerg Technical Institute
        // Aalborg University
        // November 2005
        //----------------------------------------------------------------------
        private void PrinRetMoCo(double[] SigP, double f, double[] PlasPar, double[,] Dfull, int nsigma, double[] SigP_up, out int region)
        {
            double k, m, comp, den, t1, t2, pI_II, pI_III, den_t, apex;
            double[,] D = new double[3, 3];
            double[] Rp = new double[3];
            double[] RpHat = new double[2];
            double[] Rp2 = new double[3];
            double[] Rp3 = new double[3];
            double[] N2 = new double[3];
            double[] N3 = new double[3];
            double[] NI_II = new double[3];
            double[] NI_III = new double[3];
            double[] SigPApex = new double[3];

            //--- Preliminary parameters ----------------
            k = PlasPar[0]; // Friction parameter
            m = PlasPar[2]; // Dilation parameter
            comp = PlasPar[1]; // Uniaxial compressive strength
            apex = comp / (k - 1d); // Stress coordinate of the criterions apex
            for (int i = 0; i < 3; i++)
                for (int j = 0; j < 3; j++)
                    D[i, j] = Dfull[i, j]; // Relates to normal stresses

            den = k * (D[0, 0] * m - D[0, 2]) - D[2, 0] * m + D[2, 2]; // denominator a'*D*b
            Rp[0] = (D[0, 0] * m - D[0, 2]) / den; // Rp is the scaled direction of
            Rp[1] = (D[1, 0] * m - D[1, 2]) / den; // the update direction,
            Rp[2] = (D[2, 0] * m - D[2, 2]) / den; // Rp = D*b/(a'*D*b) a, b gradient of yield surface and plastic potential, respectively.

            //--- Boundary planes -----------------------
            for (int i = 0; i < 3; i++)
                SigPApex[i] = SigP[i] - apex; // Vector from predictor stress to the apex

            // Boundary plane between regions I and II
            NI_II[0] = Rp[1] * k - Rp[2];   // NI_II = cross(Rp,R1)
            NI_II[1] = Rp[2] - Rp[0] * k; // R1 = [1 1 k], direction
            NI_II[2] = Rp[0] - Rp[1];   // vector of line 1
            pI_II = NI_II[0] * SigPApex[0] + NI_II[1] * SigPApex[1] + NI_II[2] * SigPApex[2];

            // Boundary plane between regions I and III
            NI_III[0] = Rp[1] * k - Rp[2] * k; // NI_III = cross(Rp,R2)
            NI_III[1] = Rp[2] - Rp[0] * k; // R2 = [1 k k], direction
            NI_III[2] = Rp[0] * k - Rp[1];   // vector of line 2
            pI_III = NI_III[0] * SigPApex[0] + NI_III[1] * SigPApex[1] + NI_III[2] * SigPApex[2];

            //--- t-paramters for region determination --
            // secondary surface in region II a = [0 k -1], b = [0 m -1]
            den = k * (D[1, 1] * m - D[1, 2]) - D[2, 1] * m + D[2, 2]; // denominator a'*D*b
            Rp2[0] = (D[0, 1] * m - D[0, 2]) / den; // Rp is the scaled direction of
            Rp2[1] = (D[1, 1] * m - D[1, 2]) / den; // the update direction,
            Rp2[2] = (D[2, 1] * m - D[2, 2]) / den; // Rp = D*b/(a'*D*b) a, b gradient of yield surface and plastic potential, respectively.

            N2[0] = Rp[1] * Rp2[2] - Rp[2] * Rp2[1]; // N2 = cross(Rp,Rp2)
            N2[1] = Rp[2] * Rp2[0] - Rp[0] * Rp2[2];
            N2[2] = Rp[0] * Rp2[1] - Rp[1] * Rp2[0];
            den_t = N2[0] + N2[1] + k * N2[2]; // N2'*R1, R1 = [1 1 k] Direction of line 1
                                               // t-parameter of line 1
            t1 = (N2[0] * SigPApex[0] + N2[1] * SigPApex[1] + N2[2] * SigPApex[2]) / den_t;

            // secondary surface in region III a = [k -1 0], b = [m -1 0]
            den = k * (D[0, 0] * m - D[0, 1]) - D[1, 0] * m + D[1, 1]; // denominator a'*D*b
            Rp3[0] = (D[0, 0] * m - D[0, 1]) / den; // Rp is the scaled direction of
            Rp3[1] = (D[1, 0] * m - D[1, 1]) / den; // the update direction,
            Rp3[2] = (D[2, 0] * m - D[2, 1]) / den; // Rp = D*b/(a'*D*b) a, b gradient of yield surface and plastic potential, respectively.

            N3[0] = Rp[1] * Rp3[2] - Rp[2] * Rp3[1]; // N3 = cross(Rp,Rp3)
            N3[1] = Rp[2] * Rp3[0] - Rp[0] * Rp3[2];
            N3[2] = Rp[0] * Rp3[1] - Rp[1] * Rp3[0];
            den_t = N3[0] + k * N3[1] + k * N3[2]; // N3'*R2, R2 = [1 k k] Direction of line 2
                                                   // t-parameter of line 1
            t2 = (N3[0] * (SigP[0] - apex) + N3[1] * (SigP[1] - apex) + N3[2] * (SigP[2] - apex)) / den_t;

            //--- Region determination and update -------
            if (t1 > 0d && t2 > 0d)
            {
                region = 4;
                for (int i = 0; i < 3; i++)
                    SigP_up[i] = apex;
            }
            else if (pI_II < 0d)
            {
                region = 2;
                SigP_up[0] = t1 + apex; // SigP_up = t1*R1 + apex
                SigP_up[1] = t1 + apex; // R1 = [1 1 k], direction
                SigP_up[2] = t1 * k + apex; // vector of line 1
            }
            else if (pI_III <= 0d)
            {
                region = 1;
                SigP_up[0] = SigP[0] - f * Rp[0]; // SigP_up = SigP - SiPla
                SigP_up[1] = SigP[1] - f * Rp[1]; // SiPla is the plastic corrector
                SigP_up[2] = SigP[2] - f * Rp[2]; // given by f*Rp
            }
            else
            {
                region = 3;
                SigP_up[0] = t2 + apex; // SigP_up = t2*R2 + apex
                SigP_up[1] = t2 * k + apex; // R2 = [1 k k], direction
                SigP_up[2] = t2 * k + apex; // vector of line 2
            }
        }

        //--------------------------------------------------------------------------
        // Calculates the modified elasto-plastic constitutive matrix at on a line
        // for a perfectly plastic material, It is modified in the sence that is not
        // the double singular matrix that theory prescribes. Only one type of
        // modification is possible. The calculation are carried out in principal
        // stress space.
        //
        // INPUT
        // - Name -		-type and description --
        //	mod_type	(integer,scalar) Detemines the modified type
        //				  mod_type= 0: A Doubly singular matrix, as given in the
        //							   References [1] and [2].
        //				  mod_type= 1: A single singular matrix in the plastic
        //							   strain direction, but almost singular in
        //							   any direction perpendicular to the potential
        //							   line. See Reference [3] for a detailed explanation.
        //	beta		(real,scalar) The factor by which to divide the stiffness in
        //				  the directions perpendicular to the line.
        //	SiPla		(real,array(3)) The plastic corrector principal stresses.
        //	Fline		(real,array(3)) Direction of the yield surface line in principal
        //				  stress space.
        //	Gline		(real,array(3)) Direction of the plastic potential line in principal
        //				  stress space.
        //	Dninv		(real,arrey(3,3)) Normal components of the elastic
        //				  compliance matrix.
        //	Dc			(real,array(nsigma,nsigma)) Modified elastic stiffness.
        //				  If the yield criterion is linear Dc = D.
        //	Dcinv		(real,array(nsigma,nsigma)) Inverse of the modified
        //				  elastic stiffness. If the yield criterion is linear,
        //				  Dcinv = Dinv.
        //	nsigma		(integer, scalar) Number of stress components
        //
        // OUTPUT
        //	Depc	(real,array(nsigma,nsigma) Modified Elasto-plastic constitutive
        //			  matrix on a line. If the criterion is linear Depc is the
        //			  infinitesimal version, otherwise it is the consistent. In
        //			  both cases it is expressed in principal stres space.
        //
        // REFERENCES:
        // [1] Johan Clausen, Lars Damkilde and Lars Andersen: "Efficient return
        //	  algorithms for associated plasticity with multiple yield planes", 
        //	  International Journal for Numerical Methods in Engineering, 2006,
        //	  volume 66, pages 1036-1059.
        // [2] Johan Clausen, Lars Damkilde and Lars Andersen: "An efficient return
        //	  algorithm for non-associated plasticity with linear yield criteria
        //	  in principal stress space". Computers & Structures, 2007, volume 85,
        //	  pages 1795-1807.
        // [3] Johan Clausen: Efficient Non-Linear Finite Element Implementation of
        //	  Elasto-Plasticity for Geotechnical Problems. Ph.D. Thesis. 2006.
        //	  Esbjerg Technical Institute, Aalborg University. Can be downloaded
        //	  from http://vbn.aau.dk/fbspretrieve/14058639/JCthesis.pdf
        //
        //------------------------------------------------------------------------
        // Johan Clausen
        // Department of Civil Engineering
        // Aalborg University
        // April 2008
        //------------------------------------------------------------------------
        private void FormModDepLine(int mod_type, double beta, double[] SiPla, double[] Fline, double[] Gline, double[,] Dninv, double[,] Dc,
            double[,] Dcinv, int nsigma, double[,] Depc)
        {
            double[] KoitDir = new double[3];
            double[] KoitPerDir = new double[3];
            double[,] Dper = new double[3, 3];

            Array.Clear(Depc, 0, nsigma * nsigma);
            double SiPlaDet = 0;
            for (int i = 0; i < 3; i++)
                SiPlaDet += SiPla[i] * SiPla[i];
            FormDepLinePerfect(Dc, Dcinv, Fline, Gline, nsigma, Depc);

            if (mod_type == 1 && SiPlaDet > 0d)
            {
                matmul(Dninv, SiPla, KoitDir); // the plastic strain direction
                cross(KoitDir, Gline, KoitPerDir); // Direction perpendicular to the strain direction and the plastic potential line
                double[,] DcCopy = new double[3, 3];
                double[,] DcCopyInv = new double[3, 3];
                for (int i = 0; i < 3; i++)
                    for (int j = 0; j < 3; j++)
                    {
                        DcCopy[i, j] = Dc[i, j];
                        DcCopyInv[i, j] = Dcinv[i, j];
                    }
                FormDepLinePerfect(DcCopy, DcCopyInv, KoitPerDir, KoitPerDir, 3, Dper);
                for (int i = 0; i < 3; i++)
                    for (int j = 0; j < 3; j++)
                        Depc[i, j] = Depc[i, j] + Dper[i, j] / beta;
            }
            else if (mod_type == 2 && SiPlaDet > 0d)
            {
                matmul(Dninv, SiPla, KoitDir); // The plastic strain direction in principal stress space
                double[,] DcCopy = new double[3, 3];
                double[,] DepcCopy = new double[3, 3];
                for (int i = 0; i < 3; i++)
                    for (int j = 0; j < 3; j++)
                        DcCopy[i, j] = Dc[i, j];
                FormDepPerfect(DcCopy, KoitDir, KoitDir, 3, DepcCopy);
                for (int i = 0; i < 3; i++)
                    for (int j = 0; j < 3; j++)
                        Depc[i, j] = DepcCopy[i, j];
                for (int i = 3; i < nsigma; i++)
                    for (int j = 3; j < nsigma; j++)
                        Depc[i, j] = Dc[i, j];
            }
        }

        //--------------------------------------------------------------------------
        // Calculates the modified elasto-plastic constitutive matrix at the apex
        // for a perfectly plastic material, It is modified in the sence that is not
        // the zero matrix. Two different modifications are possible, and they are
        // determined by the variable mod_type. The calculations are carried out in
        // principal stress space.
        //
        // INPUT
        // - Name -		-type and description --
        //	mod_type	(integer,scalar) Detemines the modified type
        //				  mod_type= 0: Depc is the zero matrix.
        //				  mod_type= 1: A single singular matrix in the plastic
        //							   strain direction.
        //				  mod_type= 2: A double singular matrix in the plastic
        //							   strain direction and in the direction of
        //							   the hydrostatic stress line. This should 
        //							   be used for apex points not on the hydrostatic
        //							   line (ishydro = 0)
        //	alpha		(real,scalar) The facto by which to divide the stiffness
        //	SiPla		(real,array(3)) The plastic corrector principal stresses.
        //	Dninv		(real,arrey(3,3)) Normal components of the elastic
        //				  compliance matrix.
        //	Dc			(real,array(nsigma,nsigma)) Modified elastic stiffness.
        //				  If the yield criterion is linear Dc = D.
        //	Dcinv		(real,array(nsigma,nsigma)) Inverse of the modified
        //				  elastic stiffness. If the yield criterion is linear,
        //				  Dcinv = Dinv.
        //	nsigma		(integer, scalar) Number of stress components
        //
        // OUTPUT
        //	Depc	(real,array(nsigma,nsigma) Modified Elasto-plastic constitutive
        //			  matrix on an apex. If the criterion is linear Depc is the
        //			  infinitesimal version, otherwise it is the consistent. In
        //			  both cases it is expressed in principal stres space.
        //
        //------------------------------------------------------------------------
        // Johan Clausen
        // Departmen of Civil Engineering
        // Aalborg University
        // April 2008
        //------------------------------------------------------------------------
        private void FormModDepApex(int mod_type, double alpha, double[] SiPla, double[,] Dninv, double[,] Dc, double[,] Dcinv, int nsigma, double[,] Depc)
        {
            double[] KoitDir = new double[3];
            double[] PerDir = new double[3];
            double[,] Dnorm = new double[3, 3];

            Array.Clear(Depc, 0, nsigma * nsigma);
            double SiPlaDet = 0;
            for (int i = 0; i < 3; i++)
                SiPlaDet += SiPla[i] * SiPla[i];

            if (SiPlaDet > 0d)
            {
                matmul(Dninv, SiPla, KoitDir);

                //if (mod_type == 0) then the Zero matrix is used
                if (mod_type == 1) // Single singular matrix (DepKo/fak)
                {
                    double[,] DcCopy = new double[3, 3];
                    for (int i = 0; i < 3; i++)
                        for (int j = 0; j < 3; j++)
                            DcCopy[i, j] = Dc[i, j];
                    FormDepPerfect(DcCopy, KoitDir, KoitDir, 3, Dnorm);
                    for (int i = 0; i < 3; i++)
                        for (int j = 0; j < 3; j++)
                            Depc[i, j] = Dnorm[i, j] / alpha; // Normal part
                    for (int i = 3; i < nsigma; i++)
                        for (int j = 3; j < nsigma; j++)
                            Depc[i, j] = Dc[i, j] / alpha; // Shear part
                }
                else if (mod_type == 2) // Double singular matrix (KoiteP)
                {
                    double[] Pdir = new double[] { 1d, 1d, 1d }; // Hydrostatic direction
                    cross(KoitDir, Pdir, PerDir);
                    FormDepLinePerfect(Dc, Dcinv, PerDir, PerDir, nsigma, Depc);
                    for (int i = 0; i < nsigma; i++)
                        for (int j = 0; j < nsigma; j++)
                            Depc[i, j] = Depc[i, j] / alpha; // Shear part
                }
            }
        }

        //--------------------------------------------------------------------------
        // Forms the shear part of the modification matrix T in principal stress space
        // for at perfectly plastic material. T is used when calculation the modified
        // elastic stiffness matrix Dc which is, in turn, used to form the consistent
        // constitutive matrix. The shear part of T is unaffected whether the yield
        // criterion is linear or not.
        // 
        //
        // INPUT
        //  - Name -	-type,size -	-- Description --
        //	SigP	 real,ar(3)		Principal predictor stresses in descending order
        //							  SigP = [sigP_1 sigP_2 sigP_3]
        //	SigP_up	 real,ar(3)		Updated principal stresses in descending order
        //	nshear	 int,sc			Number of shear stress components. 1 in plane problems
        //							  and 3 in full 3D
        //	s1		 int,sc			Position of largest principal in-plane stress in SigP.
        //							  Only used in plane situations
        //	s2		 int,sc			Position of smallest principal in-plane stress in SigP.
        //							  Only used in plane situations
        //
        // OUTPUT
        //	Tshear		real,ar(nshear,nshear) Shear part of the modification matrix T
        //	Tshearinv	real,ar(nshear,nshear) Inverse of Tshear
        //
        //------------------------------------------------------------------------
        // Johan Clausen
        // Department of Civil Engineering
        // Aalborg University
        // April 2008
        //------------------------------------------------------------------------
        private void TshearPrinPerfect(double[] SigP, double[] SigP_up, int nshear, int s1, int s2, double[,] Tshear, double[,] Tshearinv)
        {
            Array.Clear(Tshear, 0, nshear * nshear);
            Array.Clear(Tshearinv, 0, nshear * nshear);

            if (nshear == 1)
            {
                if (SigP_up[s1] - SigP_up[s2] > 0d)
                {
                    Tshear[0, 0] = (SigP_up[s1] - SigP_up[s2]) / (SigP[s1] - SigP[s2]);
                    Tshearinv[0, 0] = 1d / Tshear[0, 0];
                }
            }
            else if (nshear == 3)
            {
                if (SigP_up[0] - SigP_up[1] > 0d && SigP[0] - SigP[1] > 0d)
                {
                    Tshear[0, 0] = (SigP_up[0] - SigP_up[1]) / (SigP[0] - SigP[1]);
                    Tshearinv[0, 0] = 1d / Tshear[0, 0];
                }
                if (SigP_up[0] - SigP_up[2] > 0d && SigP[0] - SigP[2] > 0d)
                {
                    Tshear[1, 1] = (SigP_up[0] - SigP_up[2]) / (SigP[0] - SigP[2]);
                    Tshearinv[1, 1] = 1d / Tshear[1, 1];
                }
                if (SigP_up[1] - SigP_up[2] > 0d && SigP[1] - SigP[2] > 0d)
                {
                    Tshear[2, 2] = (SigP_up[1] - SigP_up[2]) / (SigP[1] - SigP[2]);
                    Tshearinv[2, 2] = 1d / Tshear[2, 2];
                }
            }
        }

        //--------------------------------------------------------------------------
        // Calculates the principal stress of a stress vector Sigma. Analytical 
        // expressions are used in the calculations.
        //
        // INPUT
        //  - Name - -type,size -  -- description --
        //     Sigma (real,ar(nsigma)Stress vector. Size is dependent on the type of stress
        //                   nsigma = 4 in plane situations and nsigma = 6 in 3D
        //                   The ordering of the components must be as follows:
        //                   Plane situations (plane stress, plane strain and axisymmetry):
        //                   Sigma = [sig_x sig_y sig_z tau_xy] (i.e. sig_z is out-of-plane)
        //                   General 3D:  Sigma = [sig_x sig_y sig_z tau_xy tau_xz tau_yz]
        //     nsigma  (integer,sc)  Number of stress components, i.e. the size of Sigma
        //
        // OUTPUT
        //     SigP  (real,ar(3))  Principal stresses in descending order
        //
        //---------------------------------------------------------------------------
        // Johan Clausen
        // Esbjerg Institute of Technology
        // Aalborg University
        // November 2005
        //--------------------------------------------------------------------------
        private void PrinStressAna(double[] Sigma, int nsigma, double[] SigP)
        {
            double sig_av, sig_hj, I1, J2, J3, lode, sigm, sqJ2, sin3lode;
            double[] S = new double[6];

            Array.Clear(SigP, 0, 3);
            //----- Plane situations including axisymmetry -----------------------------
            if (nsigma == 4)
            {
                sig_av = 0.5 * (Sigma[0] + Sigma[1]);
                sig_hj = Math.Sqrt(Math.Pow(0.5 * (Sigma[0] - Sigma[1]), 2) + Sigma[3] * Sigma[3]);

                SigP[0] = sig_av + sig_hj;
                SigP[1] = sig_av - sig_hj;
                SigP[2] = Sigma[2]; // The out-of-plane stress

                //-- Sorting: SigP(1) >= SigP(2) >= SigP(3) --
                if (Sigma[2] > SigP[0])
                {
                    SigP[2] = SigP[1];
                    SigP[1] = SigP[0];
                    SigP[0] = Sigma[2];
                }
                else if (Sigma[2] > SigP[1])
                {
                    SigP[2] = SigP[1];
                    SigP[1] = Sigma[2];
                }
            }
            else if (nsigma == 6)
            {
                //   ----- Invariants ----------------------------------------------------
                Invariants(Sigma, nsigma, S, out I1, out J2, out J3, out lode, out sin3lode);
                sigm = I1 / 3d; // Hydrostatic stress

                if (Math.Sqrt(J2) < 1.0e-12 * Math.Abs(I1))
                {
                    SigP[0] = Sigma[0];
                    SigP[1] = Sigma[1];
                    SigP[2] = Sigma[2];
                }
                else
                {
                    //   ----- Principal stresses -----------------------------------------
                    sqJ2 = 1.154700538379252 * Math.Sqrt(J2);
                    //				2/sqrt(3)
                    //                                   2*pi/3
                    SigP[0] = sqJ2 * Math.Sin(lode + 2.094395102393195) + sigm;
                    SigP[1] = sqJ2 * Math.Sin(lode) + sigm;
                    SigP[2] = sqJ2 * Math.Sin(lode - 2.094395102393195) + sigm;
                }
            }
        }

        //--------------------------------------------------------------------------
        // Calculates elasto-plastic constitutive matrix of a material with 
        // non-associated and perfect plasticity. If Norm = Edir the flow rule
        // is associated
        //
        // INPUT
        //     D   (real) elastic constitutive matrix. Size = nsigma x nsigma
        //     Norm  (real) Gradient (normal) of the yield surface. Size = nsigma
        //     Edir  (real) Gradient (normal) of the plastic potential. Size = nsigma
        //     nsigma  (integer) Number of stress components
        //
        // OUTPUT
        //     Dep   (real) Elasto-plastic infinitesimal constitutive matrix
        //
        //--------------------------------------------------------------------------
        // Johan Clausen
        // Esbjerg Department of Engineering
        // Aalborg University
        // July 2005
        //--------------------------------------------------------------------------
        private void FormDepPerfect(double[,] D, double[] Norm, double[] Edir, int nsigma, double[,] Dep)
        {
            double[] Num1 = new double[nsigma];
            double[] Num2 = new double[nsigma];
            double[,] Num = new double[nsigma, nsigma];

            matmul(D, Edir, Num1);
            matmul(Norm, D, Num2);
            //Num1 = matmul(D,Edir)
            //Num2 = matmul(Norm,D)
            for (int i = 0; i < nsigma; i++)
                for (int j = 0; j < nsigma; j++)
                    Num[i, j] = Num1[i] * Num2[j];

            double den = 0;
            for (int i = 0; i < nsigma; i++)
                den += Norm[i] * Num1[i];

            for (int i = 0; i < nsigma; i++)
                for (int j = 0; j < nsigma; j++)
                    Dep[i, j] = D[i, j] - Num[i, j] / den;
        }

        //--------------------------------------------------------------------------
        // Calculates double singular elasto-plastic constitutive matrix of a
        // material with non-associated and perfect plasticity. If Ra = Rb the flow
        // rule is associated. The formula is only valid in principal stress space
        //
        // INPUT
        //     D   (real) elastic constitutive matrix. Size = nsigma x nsigma
        //     Dinv  (real) Inverted elastic constitutive matrix. Size = nsigma x nsigma
        //     Ra    (real) Direction of the line on the yield surface in principal
        //           stress space. If the line is curved Ra is the tangent
        //     Rb    (real) Direction of the line on the plastic potential in principal
        //           stress space. If the line is curved Rb is the tangent
        //
        // OUTPUT
        //     Depline (real) Double singular elasto-plastic infinitesimal
        //           constitutive matrix on a line. Size = nsigma x nsigma
        //
        //------------------------------------------------------------------------
        // Johan Clausen
        // Esbjerg Department of Engineering
        // Aalborg University
        // July 2005
        //------------------------------------------------------------------------
        private void FormDepLinePerfect(double[,] D, double[,] Dinv, double[] Ra, double[] Rb, int nsigma, double[,] DepLine)
        {
            double den = 0;
            double[] Den1 = new double[3];
            double[,] Num = new double[3, 3];
            double[,] Dprin = new double[3, 3];

            for (int i = 0; i < 3; i++)
                for (int j = 0; j < 3; j++)
                    Num[i, j] = Ra[i] * Rb[j];

            double[,] DinvTemp = new double[3, 3];
            for (int i = 0; i < 3; i++)
                for (int j = 0; j < 3; j++)
                    DinvTemp[i, j] += Dinv[i, j];
            matmul(DinvTemp, Rb, Den1);
            //den = dot_product(Ra,Den1) // denomimator, Ra'*Dinv*Rb
            for (int i = 0; i < 3; i++)
                den += Ra[i] * Den1[i];

            for (int i = 0; i < 3; i++)
                for (int j = 0; j < 3; j++)
                    Dprin[i, j] = Num[i, j] / den;

            for (int i = 0; i < nsigma; i++)
                for (int j = 0; j < nsigma; j++)
                    DepLine[i, j] = 0d;
            for (int i = 0; i < 3; i++)
                for (int j = 0; j < 3; j++)
                    DepLine[i, j] = Dprin[i, j];

            if (nsigma > 3)
            {
                DepLine[3, 3] = D[3, 3]; // In plane situations nsigma == 4
                if (nsigma == 6) // General three dimensional stress state 
                {
                    DepLine[4, 4] = D[4, 4];
                    DepLine[5, 5] = D[5, 5];
                }
            }
        }

        //--------------------------------------------------------------------------
        // Calculates the principal directions of a principal stress state
        // Analytical expressions are used in the calculations.
        //
        // INPUT
        //  - Name -	-type,size -	-- description --
        //	Sigma	(real,ar(nsigma)Stress vector. Size is dependent on the type of stress
        //							  nsigma = 4 in plane situations and nsigma = 6 in 3D
        //							  The ordering of the components must be as follows:
        //							  Plane situations (plane stress, plane strain and axisymmetry):
        //							  Sigma = [sig_x sig_y sig_z tau_xy] (i.e. sig_z is out-of-plane)
        //							  General 3D:  Sigma = [sig_x sig_y sig_z tau_xy tau_xz tau_yz]
        //	nsigma	(integer,sc)	Number of stress components, i.e. the size of Sigma
        //	SigP	(real,ar(3))	Principal stresses in descending order
        //
        // OUTPUT
        //	psi		(real,ar(3,3))	Principal angle if nsigma = 4 or
        //							  vector of normalized eigenvectors (principal directions)
        //							  if nsigma = 6. When nsigma = 4 the angle is stored
        //							  in the upper left corner, i.e. psi(1,1).
        //
        //---------------------------------------------------------------------------
        // Johan Clausen
        // Esbjerg Institute of Technology
        // Aalborg University
        // November 2005
        //--------------------------------------------------------------------------
        private void PrinDirect(double[] Sigma, int nsigma, double[] SigP, double[,] psi)
        {
            const double pi = Math.PI;
            double tol1, tol2, tol, leng;
            double[] nl = new double[3];
            double[,] normal = new double[3, 3];
            int hj1, hj2, nr;

            //----- Plane situations including axisymmetry -----------------------------

            if (nsigma == 4)
            {
                double p = 0;
                if (Sigma[0] > Sigma[1] && Sigma[3] >= 0d)
                    p = 0.5 * Math.Atan(2d * Sigma[3] / (Sigma[0] - Sigma[1]));
                else if (Sigma[0] < Sigma[1] && Sigma[3] >= 0d)
                    p = 0.5 * (pi - Math.Atan(2d * Sigma[3] / (Sigma[1] - Sigma[0])));
                else if (Sigma[0] < Sigma[1] && Sigma[3] < 0d)
                    p = 0.5 * (Math.Atan(-2 * Sigma[3] / (Sigma[1] - Sigma[0])) + pi);
                else if (Sigma[0] > Sigma[1] && Sigma[3] < 0d)
                    p = 0.5 * (2d * pi - Math.Atan(-2d * Sigma[3] / (Sigma[0] - Sigma[1])));
                else if (Sigma[0] == Sigma[1] && Sigma[3] > 0d)
                    p = 0.25 * pi;
                else if (Sigma[0] == Sigma[1] && Sigma[3] < 0d)
                    p = 0.75 * pi;
                else if (Sigma[0] == Sigma[1] && Sigma[3] == 0d)
                    p = 0;

                for (int i = 0; i < 3; i++)
                    for (int j = 0; j < 3; j++)
                        psi[i, j] = p;
            }
            else if (nsigma == 6)
            {
                tol1 = Math.Abs(1e-10 * (SigP[0] - SigP[2])); // tolerance value 1
                tol2 = Math.Abs(1e-12 * (SigP[0] + SigP[1] + SigP[2])); // tolerance value 2
                tol = Math.Max(Math.Max(tol1, tol2), 1.0e-13);

                //	-- Sigma is already expressed in principal coordinates ---
                if (Math.Abs(Sigma[3]) < tol && Math.Abs(Sigma[4]) < tol && Math.Abs(Sigma[5]) < tol)
                {
                    if (Sigma[0] >= Sigma[1] && Sigma[0] >= Sigma[2])
                    {
                        normal[0, 0] = 1d;
                        if (Sigma[1] >= Sigma[2])
                        {
                            normal[1, 1] = 1d;
                            normal[2, 2] = 1d;
                        }
                        else
                        {
                            normal[2, 1] = 1d;
                            normal[1, 2] = 1d;
                        }
                    }
                    else if (Sigma[1] >= Sigma[0] && Sigma[1] >= Sigma[2])
                    {
                        normal[1, 0] = 1d;
                        if (Sigma[0] >= Sigma[2])
                        {
                            normal[0, 1] = 1d;
                            normal[2, 2] = 1d;
                        }
                        else
                        {
                            normal[2, 1] = 1d;
                            normal[0, 2] = 1d;
                        }
                    }
                    else if (Sigma[2] >= Sigma[0] && Sigma[2] >= Sigma[1])
                    {
                        normal[2, 0] = 1d;
                        if (Sigma[0] >= Sigma[1])
                        {
                            normal[0, 1] = 1d;
                            normal[1, 2] = 1d;
                        }
                        else
                        {
                            normal[1, 1] = 1d;
                            normal[0, 2] = 1d;
                        }
                    }
                }
                //	-- Three distinct eigenvalues -----
                else if (Math.Abs(SigP[0] - SigP[1]) > tol && Math.Abs(SigP[1] - SigP[2]) > tol && Math.Abs(SigP[0] - SigP[2]) > tol)
                {
                    for (int i = 0; i < 2; i++)
                    {
                        nl[0] = (Sigma[1] - SigP[i]) * (Sigma[2] - SigP[i]) - Sigma[5] * Sigma[5];
                        nl[1] = -Sigma[3] * (Sigma[2] - SigP[i]) + Sigma[5] * Sigma[4];
                        nl[2] = Sigma[3] * Sigma[5] - (Sigma[1] - SigP[i]) * Sigma[4];
                        leng = Math.Sqrt(nl[0] * nl[0] + nl[1] * nl[1] + nl[2] * nl[2]);
                        if (leng < tol)
                            normal[i, i] = 1d;
                        else
                        {
                            normal[0, i] = nl[0] / leng;
                            normal[1, i] = nl[1] / leng;
                            normal[2, i] = nl[2] / leng;
                        }
                    }
                    //			normal(:,3) = cross(Normal(:,1),Normal(:,2))
                    double[] x1, x2, x3;
                    x1 = new double[] { normal[0, 0], normal[1, 0], normal[2, 0] };
                    x2 = new double[] { normal[0, 1], normal[1, 1], normal[2, 1] };
                    x3 = new double[] { normal[0, 2], normal[1, 2], normal[2, 2] };
                    cross(x1, x2, x3);
                    normal[0, 2] = x3[0];
                    normal[1, 2] = x3[1];
                    normal[2, 2] = x3[2];
                }
                //	-- Three equal eigenvalues ----
                else if (Math.Abs(SigP[0] - SigP[1]) < tol && Math.Abs(SigP[1] - SigP[2]) < tol && Math.Abs(SigP[0] - SigP[2]) < tol)
                {
                    normal[0, 0] = 1d;
                    normal[1, 1] = 1d;
                    normal[2, 2] = 1d;
                }
                //	-- two equal eigenvalues -----
                else
                {
                    if (Math.Abs(SigP[0] - SigP[1]) <= tol)
                    {
                        hj1 = 0;
                        hj2 = 1;
                        nr = 2;
                    }
                    else if (Math.Abs(SigP[0] - SigP[2]) <= tol)
                    {
                        hj1 = 0;
                        nr = 1;
                        hj2 = 2;
                    }
                    else
                    {
                        nr = 0;
                        hj1 = 1;
                        hj2 = 2;
                    }

                    //       --- First principal direction -------------------------            
                    nl[0] = (Sigma[1] - SigP[nr]) * (Sigma[2] - SigP[nr]) - Sigma[5] * Sigma[5];
                    nl[1] = -Sigma[3] * (Sigma[2] - SigP[nr]) + Sigma[5] * Sigma[4];
                    nl[2] = Sigma[3] * Sigma[5] - (Sigma[1] - SigP[nr]) * Sigma[4];
                    leng = Math.Sqrt(nl[0] * nl[0] + nl[1] * nl[1] + nl[2] * nl[2]);
                    normal[0, nr] = nl[0] / leng;
                    normal[1, nr] = nl[1] / leng;
                    normal[2, nr] = nl[2] / leng;

                    //       --- Second principal direction ------------------------            
                    if (normal[2, nr] > 0d || normal[2, nr] < 0d)
                    {
                        nl[0] = 1d;
                        nl[1] = 1d;
                        nl[2] = -(normal[0, nr] + normal[1, nr]) / normal[2, nr];
                    }
                    else if (normal[1, nr] > 0d || normal[1, nr] < 0d)
                    {
                        nl[0] = 1d;
                        nl[2] = 1d;
                        nl[1] = -(normal[0, nr] + normal[2, nr]) / normal[1, nr];
                    }
                    else if (normal[0, nr] > 0d || normal[0, nr] < 0d)
                    {
                        nl[1] = 1d;
                        nl[2] = 1d;
                        nl[0] = -(normal[1, nr] + normal[2, nr]) / normal[0, nr];
                    }
                    leng = Math.Sqrt(nl[0] * nl[0] + nl[1] * nl[1] + nl[2] * nl[2]);
                    normal[0, hj1] = nl[0] / leng;
                    normal[1, hj1] = nl[1] / leng;
                    normal[2, hj1] = nl[2] / leng;
                    //       --- Third principal direction -------------------------
                    double[] x1 = new double[] { normal[0, nr], normal[1, nr], normal[2, nr] };
                    double[] x2 = new double[] { normal[0, hj1], normal[1, hj1], normal[2, hj1] };
                    double[] x3 = new double[3];
                    cross(x1, x2, x3);
                    normal[0, hj2] = x3[0];
                    normal[1, hj2] = x3[1];
                    normal[2, hj2] = x3[2];
                    //call cross(normal(:,nr),normal(:,hj1),normal(:,hj2))
                    //			Normal(:,hj2) = cross(normal(:,nr),Normal(:,hj1))
                }
                Array.Copy(normal, psi, 9);
            }
        }

        private void matmul(double[,] a, double[] b, double[] c)
        {
            Array.Clear(c, 0, c.Length);
            if (a.GetLength(1) == b.Length)
            {
                for (int i = 0; i < a.GetLength(0); i++)
                    for (int j = 0; j < b.Length; j++)
                        c[i] += a[i, j] * b[j];
            }
            else
                throw new ArgumentException();
        }

        private void matmul(double[] a, double[,] b, double[] c)
        {
            Array.Clear(c, 0, c.Length);
            if (b.GetLength(0) == a.Length)
            {
                for (int i = 0; i < a.Length; i++)
                    for (int j = 0; j < a.Length; j++)
                        c[i] += a[j] * b[j, i];
            }
            else
                throw new ArgumentException();
        }

        private void matmul(double[,] a, double[,] b, double[,] c)
        {
            if (a.GetLength(0) == c.GetLength(0) && b.GetLength(1) == c.GetLength(1) && a.GetLength(1) == b.GetLength(0))
            {
                for (int i = 0; i < a.GetLength(0); i++)
                {
                    for (int j = 0; j < b.GetLength(1); j++)
                    {
                        c[i, j] = 0;
                        for (int k = 0; k < a.GetLength(1); k++) // OR k<b.GetLength(0)
                            c[i, j] = c[i, j] + a[i, k] * b[k, j];
                    }
                }
            }
            else
            {
                throw new ArgumentException();
            }
        }

        //--------------------------------------------------------------------------
        // Calculates the cross product (vector product) between the two
        // three dimensional vectors A and B
        //
        // INPUT
        //     A, B  (real) Vectors with three components
        //
        // OUTPUT
        //     C (real) A vector perpendicular to both A and B
        //
        //------------------------------------------------------------------------
        // Johan Clausen
        // Esbjerg Department of Engineering
        // Aalborg University
        // July 2005
        //------------------------------------------------------------------------
        private void cross(double[] A, double[] B, double[] C)
        {
            C[0] = A[1] * B[2] - A[2] * B[1];
            C[1] = A[2] * B[0] - A[0] * B[2];
            C[2] = A[0] * B[1] - A[1] * B[0];
        }

        //--------------------------------------------------------------------
        // Calculates transformation matrix, A, for stress transformation
        // based on either an angle (plane and axisymmetry) or direction
        // cosines (general 3D). The ordering of the shear stresses
        // must be specified, as it influences how A is built
        //
        // INPUT
        //  - Name- -type,size -   -- Description --
        //     psi   (real,ar(3,3))  Matrix of direction cosines in general
        //                   three-dimensional stress state or the rotation
        //                   angle otherwise. In 3D we have 
        //                   sig_tranformed = psi'*sig*psi, where sig
        //                   is a stress tensor
        //     nsigma  (integer,sc)  Number of stress components, i.e. size of A
        //     ouplP (integer,sc)  Only used in plan situations in which it is similar to
        //                   hoop but for the principal stress vector SigP, i.e.
        //                   it is the position of the out-of-plane principal stress.
        //
        // OUTPUT
        //     A Transformation matrix, with elements ordered according to hoop and ouplP
        //         A stress is transformed according to ("'" signifies matrix transpose)
        //         Sigma_transformed = inv(A')*Sigma and conversely
        //         Sigma = A'*Sigma_transformed
        //         A strain vector, Epsilon, is tranformed according to
        //         Epsilon_transformed = A*Epsilon and conversely
        //         Epsilon = inv(A)*Epsilon_transformed
        //         See e.g. [Cook, Malkus and Plesha, 1989] for further details
        //--------------------------------------------------------------------
        // Johan Clausen
        // Esbjerg Department of Engineering
        // Aalborg University
        // July 2005
        //--------------------------------------------------------------------
        private void TransMatrix(double[,] psi, int nsigma, int ouplP, double[,] A)
        {
            double sin_psi, cos_psi, sin_psi2, cos_psi2, sin_2psi;
            double[,] A1, A2, A3, A4, Asmall, Aint;
            int[] Hj1, Hj2;
            double[,] Ext, ExtP;

            //---- Plane situations --------------------------------------------
            // Assumes that the stress components are ordered according to
            // Sigma = [sig_x sig_y sig_z tau_xy]
            //------------------------------------------------------------------ 
            A1 = new double[3, 3];
            A2 = new double[3, 3];
            A3 = new double[3, 3];
            A4 = new double[3, 3];
            Asmall = new double[3, 3];
            Aint = new double[3, 4];
            Ext = new double[3, 4];
            ExtP = new double[4, 3];

            if (nsigma == 4)
            {
                sin_psi = Math.Sin(psi[0, 0]);
                cos_psi = Math.Cos(psi[0, 0]);
                sin_psi2 = sin_psi * sin_psi;
                cos_psi2 = cos_psi * cos_psi;
                sin_2psi = Math.Sin(2d * psi[0, 0]);

                Asmall[0, 0] = cos_psi2; Asmall[0, 1] = sin_psi2;
                Asmall[1, 0] = sin_psi2; Asmall[1, 1] = cos_psi2;
                Asmall[2, 0] = -sin_2psi; Asmall[2, 1] = sin_2psi;

                Asmall[0, 2] = 0.5 * sin_2psi;
                Asmall[1, 2] = -0.5 * sin_2psi;
                Asmall[2, 2] = cos_psi2 - sin_psi2;

                if (ouplP == 2)
                {
                    ExtP[0, 0] = 1d;
                    ExtP[2, 1] = 1d;
                    ExtP[3, 2] = 1d;
                }
                else if (ouplP == 1)
                {
                    ExtP[1, 0] = 1d;
                    ExtP[2, 1] = 1d;
                    ExtP[3, 2] = 1d;
                }
                else if (ouplP == 3)
                {
                    ExtP[0, 0] = 1d;
                    ExtP[1, 1] = 1d;
                    ExtP[3, 2] = 1d;
                }

                // IMPORTANT! This order of Ext assumes that the out-of-plane stress i in position 3!!
                Ext[0, 0] = 1d;
                Ext[1, 1] = 1d;
                Ext[2, 3] = 1d;
                //Aint = matmul(Asmall,Ext)
                matmul(Asmall, Ext, Aint);
                //A = matmul(ExtP,Aint) // Transformation matrix
                matmul(ExtP, Aint, A);
                A[ouplP - 1, 2] = 1d;
            }
            //---- General three dimensional stress state ------------------------
            // Assumes that the stress components are ordered according to
            // Sigma = [sig_x sig_y sig_z tau_xy tau_xz tau_yz]
            //------------------------------------------------------------------ 
            if (nsigma == 6)
            {
                Hj1 = new int[] { 0, 2, 1 };
                Hj2 = new int[] { 1, 0, 2 };

                for (int i = 0; i < 3; i++)
                    for (int j = 0; j < 3; j++)
                        A1[i, j] = psi[i, j] * psi[i, j];

                for (int i = 0; i < 3; i++)
                    for (int j = 0; j < 3; j++)
                    {
                        A2[i, j] = 2 * psi[i, Hj1[j]] * psi[i, Hj2[j]];
                        A3[i, j] = psi[Hj1[i], j] * psi[Hj2[i], j];
                        A4[i, j] = psi[Hj1[i], Hj1[j]] * psi[Hj2[i], Hj2[j]] + psi[Hj2[i], Hj1[j]] * psi[Hj1[i], Hj2[j]];
                    }

                //A2 = 2*A2
                double[,] Atemp = new double[6, 6];
                for (int i = 0; i < 3; i++)
                    for (int j = 0; j < 3; j++)
                    {
                        Atemp[i, j] = A1[i, j];
                        Atemp[i, j + 3] = A2[i, j];
                        Atemp[i + 3, j] = A3[i, j];
                        Atemp[i + 3, j + 3] = A4[i, j];
                    }
                for (int i = 0; i < 6; i++)
                    for (int j = 0; j < 6; j++)
                        A[i, j] = Atemp[j, i];
            }
        }

        //--------------------------------------------------------------------------
        // Calculates the deviator stress vector of a stress vector Sigma, along
        // with several stress invariants.
        //
        // INPUT
        //  - Name -	-- description --
        //	Sigma	(real,ar(nsigma): Stress vector. Size is dependent on the type of stress
        //			  as given by nsigma.
        //			  nsigme = 3: It is assumed that Sigma is principal stresses,
        //						  Sigma = [SigP_1 SigP_2 SigP_3]
        //			  nsigma = 4: Plane situations, i.e. plane stress, plane strain
        //						  and axisymmetry. The ordering of the stresses must
        //						  follow Sigma = [sig_x sig_y sig_z tau_xy] , i.e.
        //						  sig_z is the out-of-plane stress.
        //			  nsigma = 6  General 3D. The stresses must be ordered according
        //						  to Sigma = [sig_x sig_y sig_z tau_xy tau_xz tau_yz]
        //	nsigma	(integer,sc): Number of stress components, i.e. the size of Sigma.
        //
        // OUTPUT
        //	S		(real, ar(nsigma)): Deviator stress vector. Same size as Sigma
        //	I1		(real,scalar): First stress invariant.
        //	J2		(real,scalar): Second deviator stress invariant.
        //	J3		(real,scalar): Third deviator stress invariant.
        //	lode	(real,scalar): Lode angle in radians. In the 
        //			  range -pi/6 <= lode <= pi/6,    (pi/6 = 30 deg)
        //	sin3lode (real,scalar): sin(3*lode)
        //
        //---------------------------------------------------------------------------
        // Johan Clausen
        // Deparment of Civil Engineering
        // Aalborg University
        // May 2008
        //--------------------------------------------------------------------------
        private void Invariants(double[] Sigma, int nsigma, double[] S, out double I1, out double J2, out double J3, out double lode, out double sin3lode)
        {
            I1 = 0;
            J2 = 0;
            J3 = 0;
            for (int i = 0; i < 3; i++)
                I1 += Sigma[i];

            if (nsigma == 3) // ! Sigma is the principal stresses
            {
                for (int i = 0; i < nsigma; i++)
                    S[i] = Sigma[i] - I1 * 0.333333333333333333; // Deviator stress vector
                J2 = 0.5 * (S[0] * S[0] + S[1] * S[1] + S[2] * S[2]);
                J3 = S[0] * S[1] * S[2];
            }
            if (nsigma == 4) // ! Plane situation, including axisymmetry
            {
                for (int i = 0; i < 3; i++)
                    S[i] = Sigma[i] - I1 * 0.333333333333333333; // Deviator stress vector
                S[3] = Sigma[3];
                J2 = 0.5 * (S[0] * S[0] + S[1] * S[1] + S[2] * S[2]) + S[3] * S[3];
                J3 = (S[0] * S[0] * S[0] + S[1] * S[1] * S[1] + S[2] * S[2] * S[2] + 3d * S[3] * S[3] * (S[0] + S[1])) * 0.3333333333333333;
            }
            if (nsigma == 6) // General three-dimensional stress state
            {
                for (int i = 0; i < 3; i++)
                    S[i] = Sigma[i] - I1 * 0.333333333333333333; // Deviator stress vector
                for (int i = 3; i < 6; i++)
                    S[i] = Sigma[i];

                J2 = 0.5 * (S[0] * S[0] + S[1] * S[1] + S[2] * S[2]) + (S[3] * S[3] + S[4] * S[4] + S[5] * S[5]);
                J3 = (S[0] * S[0] * S[0] + S[1] * S[1] * S[1] + S[2] * S[2] * S[2] + 6d * S[3] * S[4] * S[5] +
                    3d * (S[0] * (S[3] * S[3] + S[4] * S[4]) + S[1] * (S[3] * S[3] + S[5] * S[5]) + S[2] * (S[4] * S[4] + S[5] * S[5]))) * 0.3333333333333333;
            }

            if (J2 > 0d)
                sin3lode = -2.598076211353316 * J3 / Math.Pow(J2, 1.5);
            else
                sin3lode = 0d;

            if (sin3lode <= -1d)
                lode = -0.5235987755982988; // -pi/6
            else if (sin3lode >= 1d)
                lode = +0.5235987755982988; // +pi/6
            else
                lode = Math.Asin(sin3lode) / 3d;
        }
    }
}
