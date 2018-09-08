using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using ISAAR.MSolve.FEM.Interfaces;
using ISAAR.MSolve.Numerical.LinearAlgebra.Interfaces;
using System.Runtime.InteropServices;
using ISAAR.MSolve.Numerical.LinearAlgebra;
using ISAAR.MSolve.FEM.Elements.SupportiveClasses;
using ISAAR.MSolve.FEM.Embedding;
using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.Materials.Interfaces;
using ISAAR.MSolve.FEM;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.Discretization;

namespace ISAAR.MSolve.FEM.Elements
{
    public class Hexa8NL_v2 : IStructuralFiniteElement, IEmbeddedHostElement
    {
        //metavlhtes opws sto hexa8
        protected readonly static DOFType[] nodalDOFTypes = new DOFType[] { DOFType.X, DOFType.Y, DOFType.Z };
        protected readonly static DOFType[][] dofTypes = new DOFType[][] { nodalDOFTypes, nodalDOFTypes, nodalDOFTypes,
            nodalDOFTypes, nodalDOFTypes, nodalDOFTypes, nodalDOFTypes, nodalDOFTypes };
        protected readonly IContinuumMaterial3D[] materialsAtGaussPoints;
        protected IElementDOFEnumerator dofEnumerator = new GenericDOFEnumerator();
        // ews edw

        public int gp_d1_disp { get; set; } // den prepei na einai static--> shmainei idio gia ola taantikeimena afthw ths klashs
        public int gp_d2_disp { get; set; }
        public int gp_d3_disp { get; set; }
        private int nGaussPoints;

        protected Hexa8NL_v2()//consztructor apo to hexa8
        {
        }

        public Hexa8NL_v2(IContinuumMaterial3D material, int gp_d1c, int gp_d2c, int gp_d3c)
        {
            this.gp_d1_disp = gp_d1c;
            this.gp_d2_disp = gp_d2c;
            this.gp_d3_disp = gp_d3c;
            this.nGaussPoints = this.gp_d1_disp * this.gp_d2_disp * this.gp_d3_disp;
            materialsAtGaussPoints = new IContinuumMaterial3D[nGaussPoints];
			for (int i = 0; i < nGaussPoints; i++)
				materialsAtGaussPoints[i] = (IContinuumMaterial3D)material;
		}

        public int endeixiShapeFunctionAndGaussPointData = 1;
        //private double[] a_123g;
        //private double a_1g;
        //private double a_2g;
        //private double a_3g;
        //private double ksi;
        //private double heta;
        //private double zeta;
        //private int npoint;
        //private double[,] Ni;
        //private double[,] Ni_ksi;
        //private double[,] Ni_heta;
        //private double[,] Ni_zeta;
        //private double [] [,] ll1_hexa;
        //private double[][,] BL13_hexa;


        private double[][] ox_i; //den einai apo afta pou orizei o xrhsths 8 arrays twn 3 stoixeiwn
        private double[][] tx_i;
        private double[][] tu_i;
        //private double[][,] J_0b_hexa; // exoume tosa [,] osa einai kai ta gpoints
        //private double[][,] J_0_hexa;
        //private double[][,] J_0inv_hexa;
        //private double[] detJ_0; //osa kai ta gpoints
        private double[] sunt_oloklhrwmatos; // omoiws
        //private double[][,] BL11a_hexa; // exoume tosa [,] osa einai kai ta gpoints
        //private double[][,] BL12_hexa;
        //private double[][,] BL01_hexa;
        //private double[][,] BNL1_hexa;
        //private double[][,] BNL_hexa;


        //private double[,] ll2; //einai anexarthto twn GP // initialize gia to update coordinate
        private double[][] GLvec;
        private double[][] Spkvec;
        //private double[][] sunt_ol_Spkvec;
        //private double[][,] BL;
        private double[][,] Cons_disp; // upologismos apo to updatematerial


        private Tuple<double[,], double[,], double[,], double[]> GetShapeFunctionDerivativesData(int gp_d1_disp, int gp_d2_disp, int gp_d3_disp)
        {
            double[] a_123g;
            double a_1g = 0;
            double a_2g = 0;
            double a_3g = 0;
            double ksi = 0;
            double heta = 0;
            double zeta = 0;
            int npoint;
            //double[,] Ni;
            double[,] Ni_ksi;
            double[,] Ni_heta;
            double[,] Ni_zeta;

            nGaussPoints = gp_d1_disp * gp_d2_disp * gp_d3_disp;
            a_123g = new double[nGaussPoints];

            //Ni = new double[8, nGaussPoints]; // den sxetizetai me ta coh elements alla
            Ni_ksi = new double[8, nGaussPoints]; // me to prokat_disp (einai shapefunctionData)
            Ni_heta = new double[8, nGaussPoints]; // 
            Ni_zeta = new double[8, nGaussPoints];

            for (int l = 0; l < gp_d3_disp; l++)
            {
                for (int k = 0; k < gp_d2_disp; k++)
                {
                    for (int j = 0; j < gp_d1_disp; j++)
                    {
                        npoint = l * (gp_d1_disp * gp_d2_disp) + k * gp_d1_disp + j;
                        if (gp_d1_disp == 3)
                        {
                            ksi = 0.5 * (j - 1) * (j - 2) * (-0.774596669241483) + (-1) * (j) * (j - 2) * (0) + 0.5 * (j) * (j - 1) * (0.774596669241483);
                            a_1g = 0.5 * (j - 1) * (j - 2) * (0.555555555555555) + (-1) * (j) * (j - 2) * (0.888888888888888) + 0.5 * (j) * (j - 1) * (0.555555555555555);
                        }
                        if (gp_d1_disp == 2)
                        {
                            ksi = (-0.577350269189626) * (j - 1) * (-1) + (0.577350269189626) * (j) * (+1);
                            a_1g = 1;
                        }
                        if (gp_d2_disp == 3)
                        {
                            heta = 0.5 * (k - 1) * (k - 2) * (-0.774596669241483) + (-1) * (k) * (k - 2) * (0) + 0.5 * (k) * (k - 1) * (0.774596669241483);
                            a_2g = 0.5 * (k - 1) * (k - 2) * (0.555555555555555) + (-1) * (k) * (k - 2) * (0.888888888888888) + 0.5 * (k) * (k - 1) * (0.555555555555555);
                        }
                        if (gp_d2_disp == 2)
                        {
                            heta = (-0.577350269189626) * (k - 1) * (-1) + (0.577350269189626) * (k) * (+1);
                            a_2g = 1;
                        }
                        if (gp_d3_disp == 3)
                        {
                            zeta = 0.5 * (l - 1) * (l - 2) * (-0.774596669241483) + (-1) * (l) * (l - 2) * (0) + 0.5 * (l) * (l - 1) * (0.774596669241483);
                            a_3g = 0.5 * (l - 1) * (l - 2) * (0.555555555555555) + (-1) * (l) * (l - 2) * (0.888888888888888) + 0.5 * (l) * (l - 1) * (0.555555555555555);
                        }
                        if (gp_d3_disp == 2)
                        {
                            zeta = (-0.577350269189626) * (l - 1) * (-1) + (0.577350269189626) * (l) * (+1);
                            a_3g = 1;
                        }
                        a_123g[npoint] = a_1g * a_2g * a_3g;

                        //Ni[0, npoint] = 0.125 * (1 + ksi) * (1 + heta) * (1 + zeta); //N1(ksi,heta,zeta) sto sugkekrimeno GP "npoint" sto [,]
                        //Ni[1, npoint] = 0.125 * (1 - ksi) * (1 + heta) * (1 + zeta); //N2(ksi,heta,zeta) sto sugkekrimeno GP "npoint" sto [,]
                        //Ni[2, npoint] = 0.125 * (1 - ksi) * (1 - heta) * (1 + zeta);
                        //Ni[3, npoint] = 0.125 * (1 + ksi) * (1 - heta) * (1 + zeta);
                        //Ni[4, npoint] = 0.125 * (1 + ksi) * (1 + heta) * (1 - zeta);
                        //Ni[5, npoint] = 0.125 * (1 - ksi) * (1 + heta) * (1 - zeta);
                        //Ni[6, npoint] = 0.125 * (1 - ksi) * (1 - heta) * (1 - zeta);
                        //Ni[7, npoint] = 0.125 * (1 + ksi) * (1 - heta) * (1 - zeta);

                        Ni_ksi[0, npoint] = +0.125 * (1 + heta) * (1 + zeta); //N1_ksi(ksi,heta,zeta) sto sugkekrimeno GP "npoint" sto [,]
                        Ni_ksi[1, npoint] = -0.125 * (1 + heta) * (1 + zeta); //N2_ksi(ksi,heta,zeta) sto sugkekrimeno GP "npoint" sto [,]
                        Ni_ksi[2, npoint] = -0.125 * (1 - heta) * (1 + zeta);
                        Ni_ksi[3, npoint] = +0.125 * (1 - heta) * (1 + zeta);
                        Ni_ksi[4, npoint] = +0.125 * (1 + heta) * (1 - zeta);
                        Ni_ksi[5, npoint] = -0.125 * (1 + heta) * (1 - zeta);
                        Ni_ksi[6, npoint] = -0.125 * (1 - heta) * (1 - zeta);
                        Ni_ksi[7, npoint] = +0.125 * (1 - heta) * (1 - zeta);

                        Ni_heta[0, npoint] = 0.125 * (1 + ksi) * (+1) * (1 + zeta); //N1_heta(ksi,heta,zeta) sto sugkekrimeno GP "npoint" sto [,]
                        Ni_heta[1, npoint] = 0.125 * (1 - ksi) * (+1) * (1 + zeta); //N2_heta(ksi,heta,zeta) sto sugkekrimeno GP "npoint" sto [,]
                        Ni_heta[2, npoint] = 0.125 * (1 - ksi) * (-1) * (1 + zeta);
                        Ni_heta[3, npoint] = 0.125 * (1 + ksi) * (-1) * (1 + zeta);
                        Ni_heta[4, npoint] = 0.125 * (1 + ksi) * (+1) * (1 - zeta);
                        Ni_heta[5, npoint] = 0.125 * (1 - ksi) * (+1) * (1 - zeta);
                        Ni_heta[6, npoint] = 0.125 * (1 - ksi) * (-1) * (1 - zeta);
                        Ni_heta[7, npoint] = 0.125 * (1 + ksi) * (-1) * (1 - zeta);

                        Ni_zeta[0, npoint] = 0.125 * (1 + ksi) * (1 + heta) * (+1); //N1_zeta(ksi,heta,zeta) sto sugkekrimeno GP "npoint" sto [,]
                        Ni_zeta[1, npoint] = 0.125 * (1 - ksi) * (1 + heta) * (+1); //N2_zeta(ksi,heta,zeta) sto sugkekrimeno GP "npoint" sto [,]
                        Ni_zeta[2, npoint] = 0.125 * (1 - ksi) * (1 - heta) * (+1);
                        Ni_zeta[3, npoint] = 0.125 * (1 + ksi) * (1 - heta) * (+1);
                        Ni_zeta[4, npoint] = 0.125 * (1 + ksi) * (1 + heta) * (-1);
                        Ni_zeta[5, npoint] = 0.125 * (1 - ksi) * (1 + heta) * (-1);
                        Ni_zeta[6, npoint] = 0.125 * (1 - ksi) * (1 - heta) * (-1);
                        Ni_zeta[7, npoint] = 0.125 * (1 + ksi) * (1 - heta) * (-1);
                    }
                }
            }

            Tuple<double[,], double[,], double[,], double[]> shapeFunctionDerivativesData = new Tuple<double[,], double[,], double[,], double[]>(Ni_ksi, Ni_heta, Ni_zeta, a_123g);
            return shapeFunctionDerivativesData;

        }

        private double[][,] Getll1Hexa(double[,] Ni_ksi, double[,] Ni_heta, double[,] Ni_zeta, int gp_d1_disp, int gp_d2_disp, int gp_d3_disp)
        {
            int nGaussPoints = gp_d1_disp * gp_d2_disp * gp_d3_disp;
            int npoint;

            double[][,] ll1_hexa;
            ll1_hexa = new double[nGaussPoints][,];
            for (int l = 0; l < gp_d3_disp; l++)
            {
                for (int k = 0; k < gp_d2_disp; k++)
                {
                    for (int j = 0; j < gp_d1_disp; j++)
                    {
                        npoint = l * (gp_d1_disp * gp_d2_disp) + k * gp_d1_disp + j;

                        ll1_hexa[npoint] = new double[3, 8];
                        for (int m = 0; m < 8; m++)
                        {
                            ll1_hexa[npoint][0, m] = Ni_ksi[m, npoint];
                            ll1_hexa[npoint][1, m] = Ni_heta[m, npoint];
                            ll1_hexa[npoint][2, m] = Ni_zeta[m, npoint];
                        }

                    }
                }
            }
            return ll1_hexa;
        }

        private double[][,] GetBL13Hexa(double[][,] ll1_hexa, int gp_d1_disp, int gp_d2_disp, int gp_d3_disp)
        {
            double[][,] BL13_hexa;
            int nGaussPoints = gp_d1_disp * gp_d2_disp * gp_d3_disp;
            int npoint;
            BL13_hexa = new double[nGaussPoints][,];
            for (int l = 0; l < gp_d3_disp; l++)
            {
                for (int k = 0; k < gp_d2_disp; k++)
                {
                    for (int j = 0; j < gp_d1_disp; j++)
                    {
                        npoint = l * (gp_d1_disp * gp_d2_disp) + k * gp_d1_disp + j;

                        BL13_hexa[npoint] = new double[9, 24];
                        for (int m = 0; m < 9; m++)
                        {
                            for (int n = 0; n < 24; n++)
                            { BL13_hexa[npoint][m, n] = 0; }
                        }// TODO: den xreiazetai mhdenismos edw 

                        for (int m = 0; m < 8; m++)
                        {
                            for (int n = 0; n < 3; n++)
                            {
                                BL13_hexa[npoint][n, 3 * m + 0] = ll1_hexa[npoint][n, m];
                                BL13_hexa[npoint][n + 3, 3 * m + 1] = ll1_hexa[npoint][n, m];
                                BL13_hexa[npoint][n + 6, 3 * m + 2] = ll1_hexa[npoint][n, m];
                            }
                        }
                    }
                }

            }
            return BL13_hexa;
        }

        private Tuple<double[][,], double[]> GetJ_0invHexaAndDetJ_0(double[][,] ll1_hexa, IElement element)
        {
            double[][,] J_0b_hexa; // exoume tosa [,] osa einai kai ta gpoints
            double[][,] J_0_hexa;
            double[][,] J_0inv_hexa;
            double[] detJ_0; //osa kai ta gpoints

            double[][] ox_i;
            ox_i = new double[8][];
            for (int j = 0; j < 8; j++)
            {
                ox_i[j] = new double[] { element.INodes[j].X, element.INodes[j].Y, element.INodes[j].Z, };
            }
            J_0b_hexa = new double[nGaussPoints][,];
            J_0_hexa = new double[nGaussPoints][,];
            J_0inv_hexa = new double[nGaussPoints][,];
            detJ_0 = new double[nGaussPoints];

            for (int gpoint = 0; gpoint < nGaussPoints; gpoint++)
            {
                // initialize diastaseis twn mhtrwwn kai meta gemisma keliwn (olwn h mono oswn mporoume sthn arxh)
                J_0b_hexa[gpoint] = new double[8, 3];
                J_0_hexa[gpoint] = new double[3, 3];
                J_0inv_hexa[gpoint] = new double[3, 3];

                //
                for (int m = 0; m < 8; m++)
                {
                    for (int n = 0; n < 3; n++)
                    {
                        J_0b_hexa[gpoint][m, n] = ox_i[m][n];
                    }
                }

                //
                for (int m = 0; m < 3; m++)
                {
                    for (int n = 0; n < 3; n++)
                    {
                        J_0_hexa[gpoint][m, n] = 0;
                        for (int p = 0; p < 8; p++)
                        {
                            J_0_hexa[gpoint][m, n] += ll1_hexa[gpoint][m, p] * J_0b_hexa[gpoint][p, n];
                        }
                    }
                }

                //
                double det1 = J_0_hexa[gpoint][0, 0] *
                         ((J_0_hexa[gpoint][1, 1] * J_0_hexa[gpoint][2, 2]) - (J_0_hexa[gpoint][2, 1] * J_0_hexa[gpoint][1, 2]));
                double det2 = J_0_hexa[gpoint][0, 1] *
                              ((J_0_hexa[gpoint][1, 0] * J_0_hexa[gpoint][2, 2]) - (J_0_hexa[gpoint][2, 0] * J_0_hexa[gpoint][1, 2]));
                double det3 = J_0_hexa[gpoint][0, 2] *
                              ((J_0_hexa[gpoint][1, 0] * J_0_hexa[gpoint][2, 1]) - (J_0_hexa[gpoint][2, 0] * J_0_hexa[gpoint][1, 1]));
                double jacobianDeterminant = det1 - det2 + det3;
                if (jacobianDeterminant < 0)
                {
                    throw new InvalidOperationException("The Jacobian Determinant is negative.");
                }
                detJ_0[gpoint] = jacobianDeterminant;

                //
                J_0inv_hexa[gpoint][0, 0] = ((J_0_hexa[gpoint][1, 1] * J_0_hexa[gpoint][2, 2]) - (J_0_hexa[gpoint][2, 1] * J_0_hexa[gpoint][1, 2])) *
                                    (1 / detJ_0[gpoint]);
                J_0inv_hexa[gpoint][0, 1] = ((J_0_hexa[gpoint][2, 1] * J_0_hexa[gpoint][0, 2]) - (J_0_hexa[gpoint][0, 1] * J_0_hexa[gpoint][2, 2])) *
                                        (1 / detJ_0[gpoint]);
                J_0inv_hexa[gpoint][0, 2] = ((J_0_hexa[gpoint][0, 1] * J_0_hexa[gpoint][1, 2]) - (J_0_hexa[gpoint][1, 1] * J_0_hexa[gpoint][0, 2])) *
                                        (1 / detJ_0[gpoint]);
                J_0inv_hexa[gpoint][1, 0] = ((J_0_hexa[gpoint][2, 0] * J_0_hexa[gpoint][1, 2]) - (J_0_hexa[gpoint][1, 0] * J_0_hexa[gpoint][2, 2])) *
                                        (1 / detJ_0[gpoint]);
                J_0inv_hexa[gpoint][1, 1] = ((J_0_hexa[gpoint][0, 0] * J_0_hexa[gpoint][2, 2]) - (J_0_hexa[gpoint][2, 0] * J_0_hexa[gpoint][0, 2])) *
                                        (1 / detJ_0[gpoint]);
                J_0inv_hexa[gpoint][1, 2] = ((J_0_hexa[gpoint][1, 0] * J_0_hexa[gpoint][0, 2]) - (J_0_hexa[gpoint][0, 0] * J_0_hexa[gpoint][1, 2])) *
                                        (1 / detJ_0[gpoint]);
                J_0inv_hexa[gpoint][2, 0] = ((J_0_hexa[gpoint][1, 0] * J_0_hexa[gpoint][2, 1]) - (J_0_hexa[gpoint][2, 0] * J_0_hexa[gpoint][1, 1])) *
                                        (1 / detJ_0[gpoint]);
                J_0inv_hexa[gpoint][2, 1] = ((J_0_hexa[gpoint][2, 0] * J_0_hexa[gpoint][0, 1]) - (J_0_hexa[gpoint][2, 1] * J_0_hexa[gpoint][0, 0])) *
                                        (1 / detJ_0[gpoint]);
                J_0inv_hexa[gpoint][2, 2] = ((J_0_hexa[gpoint][0, 0] * J_0_hexa[gpoint][1, 1]) - (J_0_hexa[gpoint][1, 0] * J_0_hexa[gpoint][0, 1])) *
                                        (1 / detJ_0[gpoint]);
            }


            Tuple<double[][,], double[]> J_0inv_hexaAndDetJ_0 = new Tuple<double[][,], double[]>(J_0inv_hexa, detJ_0);
            return J_0inv_hexaAndDetJ_0;

        }

        private double[][,] GetBL11a_hexa(double[][,] J_0inv_hexa)
        {
            double[][,] BL11a_hexa = new double[nGaussPoints][,];
            for (int gpoint = 0; gpoint < nGaussPoints; gpoint++)
            {
                BL11a_hexa[gpoint] = new double[6, 9];
                for (int m = 0; m < 3; m++) // upologismos triwn prwtwn grammwn
                {
                    for (int n = 0; n < 3; n++)
                    {
                        BL11a_hexa[gpoint][m, 3 * m + n] = J_0inv_hexa[gpoint][m, n];
                    }
                }
                for (int n = 0; n < 3; n++)
                {
                    BL11a_hexa[gpoint][3, n] = J_0inv_hexa[gpoint][1, n]; // upologismos 4hs gramms
                    BL11a_hexa[gpoint][3, 3 + n] = J_0inv_hexa[gpoint][0, n];
                    BL11a_hexa[gpoint][4, 3 + n] = J_0inv_hexa[gpoint][2, n]; // upologismos 5hs gramms
                    BL11a_hexa[gpoint][4, 6 + n] = J_0inv_hexa[gpoint][1, n];
                    BL11a_hexa[gpoint][5, 0 + n] = J_0inv_hexa[gpoint][2, n]; // upologismos 6hs gramms
                    BL11a_hexa[gpoint][5, 6 + n] = J_0inv_hexa[gpoint][0, n];
                }
            }

            return BL11a_hexa;
        }

        private double[][,] GetBL12_hexa(double[][,] J_0inv_hexa)
        {
            double[][,] BL12_hexa = new double[nGaussPoints][,];
            for (int gpoint = 0; gpoint < nGaussPoints; gpoint++)
            {
                BL12_hexa[gpoint] = new double[9, 9];
                for (int m = 0; m < 3; m++) // upologismos triwn prwtwn grammwn
                {
                    for (int n = 0; n < 3; n++)
                    {
                        BL12_hexa[gpoint][m, 3 * m + n] = J_0inv_hexa[gpoint][0, n];
                    }
                }
                for (int m = 0; m < 3; m++) // upologismos grammwn 4 ews 6
                {
                    for (int n = 0; n < 3; n++)
                    {
                        BL12_hexa[gpoint][3 + m, 3 * m + n] = J_0inv_hexa[gpoint][1, n];
                    }
                }
                for (int m = 0; m < 3; m++) // upologismos grammwn 7 ews 8
                {
                    for (int n = 0; n < 3; n++)
                    {
                        BL12_hexa[gpoint][6 + m, 3 * m + n] = J_0inv_hexa[gpoint][2, n];
                    }
                }

            }

            return BL12_hexa;
        }

        private double[][,] GetBL01_hexa(double[][,] J_0inv_hexa)
        {
            double[][,] BL01_hexa = new double[nGaussPoints][,];
            for (int gpoint = 0; gpoint < nGaussPoints; gpoint++)
            {
                BL01_hexa[gpoint] = new double[6, 9];
                for (int m = 0; m < 3; m++) // upologismos triwn prwtwn grammwn
                {
                    for (int n = 0; n < 3; n++)
                    {
                        BL01_hexa[gpoint][m, 3 * m + n] = J_0inv_hexa[gpoint][m, n];
                    }
                }
                for (int n = 0; n < 3; n++)
                {
                    BL01_hexa[gpoint][3, n] = J_0inv_hexa[gpoint][1, n]; // upologismos 4hs gramms
                    BL01_hexa[gpoint][3, 3 + n] = J_0inv_hexa[gpoint][0, n];
                    BL01_hexa[gpoint][4, 3 + n] = J_0inv_hexa[gpoint][2, n]; // upologismos 5hs gramms
                    BL01_hexa[gpoint][4, 6 + n] = J_0inv_hexa[gpoint][1, n];
                    BL01_hexa[gpoint][5, 0 + n] = J_0inv_hexa[gpoint][2, n]; // upologismos 6hs gramms
                    BL01_hexa[gpoint][5, 6 + n] = J_0inv_hexa[gpoint][0, n];
                }
            }
            return BL01_hexa;
        }

        private double[][,] GetBNL1_hexa(double[][,] J_0inv_hexa)
        {
            double[][,] BNL1_hexa = new double[nGaussPoints][,];
            for (int gpoint = 0; gpoint < nGaussPoints; gpoint++)
            {
                BNL1_hexa[gpoint] = new double[9, 9];
                for (int m = 0; m < 3; m++)
                {
                    for (int n = 0; n < 3; n++)
                    {
                        for (int p = 0; p < 3; p++)
                        {
                            BNL1_hexa[gpoint][3 * m + n, 3 * m + p] = J_0inv_hexa[gpoint][n, p];
                        }
                    }
                }
            }
            return BNL1_hexa;
        }



        private void CalculateInitialConfigurationData(IElement element)
        {
            double[] a_123g;


            //double[,] Ni;
            double[,] Ni_ksi;
            double[,] Ni_heta;
            double[,] Ni_zeta;

            Tuple<double[,], double[,], double[,], double[]> shapeFunctionDerivativesData;
            shapeFunctionDerivativesData = GetShapeFunctionDerivativesData(gp_d1_disp, gp_d2_disp, gp_d3_disp);
            Ni_ksi = shapeFunctionDerivativesData.Item1;
            Ni_heta = shapeFunctionDerivativesData.Item2;
            Ni_zeta = shapeFunctionDerivativesData.Item3;
            a_123g = shapeFunctionDerivativesData.Item4;

            double[][,] ll1_hexa;
            ll1_hexa = Getll1Hexa(Ni_ksi, Ni_heta, Ni_zeta, gp_d1_disp, gp_d2_disp, gp_d3_disp);
            double[][,] BL13_hexa;
            BL13_hexa = GetBL13Hexa(ll1_hexa, gp_d1_disp, gp_d2_disp, gp_d3_disp);

            double[][,] BNL1_hexa;

            //PrintUtilities.WriteToFile(Ni, @"C:\Users\turbo-x\Desktop\cohesive_check_MSOLVE_2\paradeigma_apo_arxika_swsta_embeded_shell_gia_check_tou_rve_embedding_sto_MSolve\elegxos_alalgwn_fe2_tax_me1_arxiko_chol_dixws_me1_OneElementRVECheckExample\Ni_data_8epi27.txt");
            endeixiShapeFunctionAndGaussPointData = 2;

            ox_i = new double[8][];
            tx_i = new double[8][];
            tu_i = new double[8][]; // apla initialized edw kai tpt allo

            Tuple<double[][,], double[]> J_0inv_hexaAndDetJ_0;
            J_0inv_hexaAndDetJ_0 = GetJ_0invHexaAndDetJ_0(ll1_hexa, element);
            double[][,] J_0inv_hexa;
            J_0inv_hexa = J_0inv_hexaAndDetJ_0.Item1;
            double[] detJ_0; //osa kai ta gpoints
            detJ_0 = J_0inv_hexaAndDetJ_0.Item2;

            sunt_oloklhrwmatos = new double[nGaussPoints];

            //BL11a_hexa = GetBL11a_hexa(J_0inv_hexa);
            //BL12_hexa = GetBL12_hexa(J_0inv_hexa);
            //BL01_hexa = GetBL01_hexa(J_0inv_hexa);
            BNL1_hexa = GetBNL1_hexa(J_0inv_hexa);

            //BNL_hexa = new double[nGaussPoints][,];

            for (int j = 0; j < 8; j++)
            {
                ox_i[j] = new double[] { element.INodes[j].X, element.INodes[j].Y, element.INodes[j].Z, };
                tx_i[j] = new double[] { element.INodes[j].X, element.INodes[j].Y, element.INodes[j].Z, };
                //tu_i[j] = new double[] { 0, 0, 0 }; den ananewnontai se afth th methodo ta mhtrwa pou periexoun tu_i
            }

            for (int gpoint = 0; gpoint < nGaussPoints; gpoint++)
            {
                // initialize diastaseis twn mhtrwwn kai meta gemisma keliwn (olwn h mono oswn mporoume sthn arxh)
                //BNL_hexa[gpoint] = new double[9, 24];

                sunt_oloklhrwmatos[gpoint] = detJ_0[gpoint] * a_123g[gpoint];

                //
                //for (int m = 0; m < 9; m++)
                //{
                //    for (int n = 0; n < 24; n++)
                //    {
                //        BNL_hexa[gpoint][m, n] = 0;
                //        for (int p = 0; p < 9; p++)
                //        {
                //            BNL_hexa[gpoint][m, n] += BNL1_hexa[gpoint][m, p] * BL13_hexa[gpoint][p, n];
                //        }
                //    }
                //}
            }



            tu_i = new double[8][];
            tx_i = new double[8][];
            //ll2 = new double[8, 3];
            GLvec = new double[nGaussPoints][];
            Spkvec = new double[nGaussPoints][];
            for (int gpoint = 0; gpoint < nGaussPoints; gpoint++)
            {
                GLvec[gpoint] = new double[6];
                Spkvec[gpoint] = new double[6];
            }
            for (int k = 0; k < 8; k++)
            {
                tu_i[k] = new double[3];
                tx_i[k] = new double[3];
            }

            //sunt_ol_Spkvec = new double[nGaussPoints][];
            //BL = new double[nGaussPoints][,];
            for (int gpoint = 0; gpoint < nGaussPoints; gpoint++)
            {
                //sunt_ol_Spkvec[gpoint] = new double[6];
                //BL[gpoint] = new double[6, 24];

            }
            Cons_disp = new double[nGaussPoints][,]; // upologismos apo to updatematerial
            for (int gpoint = 0; gpoint < nGaussPoints; gpoint++)
            {
                Cons_disp[gpoint] = new double[6, 6];
            }

        }




        private void InitializeMatrices()
        {
            // initialize gia to update coordinate
            tu_i = new double[8][];
            tx_i = new double[8][];
            //ll2 = new double[8, 3];
            //J_1b = new double[8, 3];
            //J_1 = new double[nGaussPoints][,];
            //BL11b= new double[nGaussPoints][,];
            //DGtr = new double[nGaussPoints][,];
            //GL = new double[nGaussPoints][,];
            GLvec = new double[nGaussPoints][];
            Spkvec = new double[nGaussPoints][];
            for (int npoint = 0; npoint < nGaussPoints; npoint++)
            {
                //BL11b[npoint] = new double[9,9];   
                //J_1[npoint] = new double[3, 3];
                //DGtr[npoint] = new double[3, 3];
                //GL[npoint] = new double[3, 3];
                GLvec[npoint] = new double[6];
                Spkvec[npoint] = new double[6];
                // epiprosthtws mhdenizoueme gia efkolia edw osa apo ta parapanw the tha gemisoun plhrws 
                //for (int j = 0; j < 9; j++)
                //{
                //    for (int k = 0; k < 9; k++)
                //    {
                //        BL11b[npoint][j, k] = 0;
                //    }
                //}
            }
            for (int k = 0; k < 8; k++)
            {
                tu_i[k] = new double[3];
                tx_i[k] = new double[3];
            }

            // initialize gia to update forces
            //l_perisp = new double[3, 3];
            //sunt_ol_Spkvec = new double[nGaussPoints][];
            //BL11 = new double[nGaussPoints][,];
            //BL1112sun01_hexa = new double[nGaussPoints][,];
            //BL = new double[nGaussPoints][,];
            //fxk1 = new double[nGaussPoints+1][];
            for (int npoint = 0; npoint < nGaussPoints; npoint++)
            {
                //sunt_ol_Spkvec[npoint] = new double[6];
                //BL11[npoint] = new double[6, 9];
                //BL1112sun01_hexa[npoint] = new double[6, 9];
                //BL[npoint] = new double[6, 24];

            }
            //for (int npoint = 0; npoint < nGaussPoints+1; npoint++)
            //{
            //    fxk1[npoint] = new double[24];
            //}

            // initialize gia to updateKmatrices
            //sunt_ol_Spk = new double[nGaussPoints][,];
            //sunt_ol_SPK_epi_BNL_hexa = new double[9,24];
            //kl_ = new double[nGaussPoints + 1][,]; 
            //knl_ = new double[nGaussPoints + 1][,]; 
            //k_stoixeiou = new double[24,24]; 
            Cons_disp = new double[nGaussPoints][,]; // upologismos apo to updatematerial
            //sunt_ol_cons_disp = new double [6,6];
            //sunt_ol_cons_disp_epi_BL = new double[6, 24];
            for (int npoint = 0; npoint < nGaussPoints; npoint++)
            {
                //sunt_ol_Spk[npoint] = new double[3,3];
                Cons_disp[npoint] = new double[6, 6];
            }
            //for (int npoint = 0; npoint < nGaussPoints + 1; npoint++)
            //{
            //    kl_[npoint] = new double[24,24];
            //    knl_[npoint] = new double[24, 24];
            //}


        }



        private void UpdateCoordinateData(double[] localdisplacements, IElement element) // sto shell8disp sto calculate forces kaleitai me this.UpdateCoordinateData(localTotalDisplacements);
        {
            //YPOLOGISMOS EDW KAI TOU ll1_hexa pou de tha karatietai pia kai olwn 
            double[,] Ni_ksi;
            double[,] Ni_heta;
            double[,] Ni_zeta;

            Tuple<double[,], double[,], double[,], double[]> shapeFunctionDerivativesData;
            shapeFunctionDerivativesData = GetShapeFunctionDerivativesData(gp_d1_disp, gp_d2_disp, gp_d3_disp);
            Ni_ksi = shapeFunctionDerivativesData.Item1;
            Ni_heta = shapeFunctionDerivativesData.Item2;
            Ni_zeta = shapeFunctionDerivativesData.Item3;
            double[][,] ll1_hexa;
            ll1_hexa = Getll1Hexa(Ni_ksi, Ni_heta, Ni_zeta, gp_d1_disp, gp_d2_disp, gp_d3_disp);
            Tuple<double[][,], double[]> J_0inv_hexaAndDetJ_0;
            J_0inv_hexaAndDetJ_0 = GetJ_0invHexaAndDetJ_0(ll1_hexa, element);
            double[][,] J_0inv_hexa;
            J_0inv_hexa = J_0inv_hexaAndDetJ_0.Item1;
            //YPOLOGISMOS EDW KAI TOU ll1_hexa pou de tha karatietai pia kai olwn kai tou J_0inv


            double[,] J_1b = new double[8, 3];
            double[][,] J_1 = new double[nGaussPoints][,];
            double[][,] DGtr = new double[nGaussPoints][,];
            double[][,] GL = new double[nGaussPoints][,];
            for (int npoint = 0; npoint < nGaussPoints; npoint++)
            {
                J_1[npoint] = new double[3, 3];
                DGtr[npoint] = new double[3, 3];
                GL[npoint] = new double[3, 3];
            }

            for (int j = 0; j < 8; j++)
            {
                for (int k = 0; k < 3; k++)
                {
                    tu_i[j][k] = localdisplacements[3 * j + k];
                    tx_i[j][k] = ox_i[j][k] + tu_i[j][k];
                }
            }

            //
            for (int m = 0; m < 8; m++)
            {
                for (int n = 0; n < 3; n++)
                {
                    //ll2[m, n] = tu_i[m][n];
                    J_1b[m, n] = tx_i[m][n];
                }
            }

            // //
            for (int npoint = 0; npoint < nGaussPoints; npoint++)
            {

                //
                for (int m = 0; m < 3; m++)
                {
                    for (int n = 0; n < 3; n++)
                    {
                        J_1[npoint][m, n] = 0;
                        for (int p = 0; p < 8; p++)
                        {
                            J_1[npoint][m, n] += ll1_hexa[npoint][m, p] * J_1b[p, n];
                        }
                    }
                }


                //
                //for (int m = 0; m < 3; m++)
                //{
                //    for (int n = 0; n < 3; n++)
                //    {
                //        for (int p = 0; p < 3; p++)
                //        {
                //            BL11b[npoint][3 * m + n, 3 * m + p] = l_perisp[n, p]; 
                //        }
                //    }
                //}

                //
                for (int m = 0; m < 3; m++)
                {
                    for (int n = 0; n < 3; n++)
                    {
                        DGtr[npoint][m, n] = 0;
                        for (int p = 0; p < 3; p++)
                        {
                            DGtr[npoint][m, n] += J_0inv_hexa[npoint][m, p] * J_1[npoint][p, n];
                        }
                    }
                }

                //
                for (int m = 0; m < 3; m++)
                {
                    for (int n = 0; n < 3; n++)
                    {
                        GL[npoint][m, n] = 0;
                        for (int p = 0; p < 3; p++)
                        {
                            GL[npoint][m, n] += DGtr[npoint][m, p] * DGtr[npoint][n, p];
                        }
                    }
                }
                for (int m = 0; m < 3; m++)
                {
                    GL[npoint][m, m] += -1;
                }
                for (int m = 0; m < 3; m++)
                {
                    for (int n = 0; n < 3; n++)
                    {
                        GL[npoint][m, n] = 0.5 * GL[npoint][m, n];
                    }
                }

                //
                for (int m = 0; m < 3; m++)
                {
                    GLvec[npoint][m] = GL[npoint][m, m];
                }
                GLvec[npoint][3] = 2 * GL[npoint][0, 1];
                GLvec[npoint][4] = 2 * GL[npoint][1, 2];
                GLvec[npoint][5] = 2 * GL[npoint][2, 0];
            }

        }


        // apo uliko tha einai gnwsto to SpkVec[npoint][1:6] gia ola ta npoints
        // me vash afto programmatizontai oi forces


        private double[] UpdateForces(Element element)
        {
            // upologismos entos forces olwn twn apaitoumenwn mhtrwwn
            double[,] ll2 = new double[8, 3];
            for (int m = 0; m < 8; m++)
            {
                for (int n = 0; n < 3; n++)
                {
                    ll2[m, n] = tu_i[m][n];
                }
            }
            // upologismos ennoeitai kai twn mhtrwwn apo tous arxikous upologismous pou xreiazontai edw (ola ews ll1_hexa)
            double[,] Ni_ksi;
            double[,] Ni_heta;
            double[,] Ni_zeta;
            Tuple<double[,], double[,], double[,], double[]> shapeFunctionDerivativesData;
            shapeFunctionDerivativesData = GetShapeFunctionDerivativesData(gp_d1_disp, gp_d2_disp, gp_d3_disp);
            Ni_ksi = shapeFunctionDerivativesData.Item1;
            Ni_heta = shapeFunctionDerivativesData.Item2;
            Ni_zeta = shapeFunctionDerivativesData.Item3;
            double[][,] ll1_hexa;
            ll1_hexa = Getll1Hexa(Ni_ksi, Ni_heta, Ni_zeta, gp_d1_disp, gp_d2_disp, gp_d3_disp);
            //11a 12 13 kai 01 epishs xreiazontai opote kai to J_0inv_hexa pou afta theloun 
            Tuple<double[][,], double[]> J_0inv_hexaAndDetJ_0;
            J_0inv_hexaAndDetJ_0 = GetJ_0invHexaAndDetJ_0(ll1_hexa, element);
            double[][,] J_0inv_hexa;
            J_0inv_hexa = J_0inv_hexaAndDetJ_0.Item1;
            double[][,] BL13_hexa;
            BL13_hexa = GetBL13Hexa(ll1_hexa, gp_d1_disp, gp_d2_disp, gp_d3_disp);
            double[][,] BL11a_hexa; // exoume tosa [,] osa einai kai ta gpoints
            double[][,] BL12_hexa;
            double[][,] BL01_hexa;
            BL11a_hexa = GetBL11a_hexa(J_0inv_hexa);
            BL12_hexa = GetBL12_hexa(J_0inv_hexa);
            BL01_hexa = GetBL01_hexa(J_0inv_hexa);

            //INITIALIZE EDW TWN mhtrwwn pou den tha apothikevontai pia
            double[][] sunt_ol_Spkvec = new double[nGaussPoints][];
            double[][,] BL = new double[nGaussPoints][,];
            for (int gpoint = 0; gpoint < nGaussPoints; gpoint++)
            {
                sunt_ol_Spkvec[gpoint] = new double[6];
                BL[gpoint] = new double[6, 24];

            }



            double[][] fxk1 = new double[nGaussPoints + 1][];
            for (int npoint = 0; npoint < nGaussPoints + 1; npoint++)
            {
                fxk1[npoint] = new double[24];
            }

            double[][,] BL11 = new double[nGaussPoints][,];
            double[][,] BL1112sun01_hexa = new double[nGaussPoints][,];
            for (int npoint = 0; npoint < nGaussPoints; npoint++)
            {
                BL11[npoint] = new double[6, 9];
                BL1112sun01_hexa[npoint] = new double[6, 9];
            }



            for (int npoint = 0; npoint < nGaussPoints; npoint++)
            {

                //
                for (int m = 0; m < 6; m++)
                {
                    sunt_ol_Spkvec[npoint][m] = sunt_oloklhrwmatos[npoint] * Spkvec[npoint][m];
                }

                //
                double[,] l_perisp = new double[3, 3];
                for (int m = 0; m < 3; m++)
                {
                    for (int n = 0; n < 3; n++)
                    {
                        //l_perisp[m, n] = 0; //sp1
                        for (int p = 0; p < 8; p++)
                        {
                            l_perisp[m, n] += ll1_hexa[npoint][m, p] * ll2[p, n];
                        }
                    }
                }

                //
                //for (int m = 0; m < 6; m++)
                //{
                //    for (int n = 0; n < 9; n++)
                //    {
                //        BL11[npoint][m,n] = 0;
                //    }
                //} sp1
                for (int m = 0; m < 6; m++)
                {
                    for (int n = 0; n < 3; n++)
                    {
                        for (int p = 0; p < 3; p++)
                        {
                            BL11[npoint][m, n] += BL11a_hexa[npoint][m, p] * l_perisp[p, n];
                            BL11[npoint][m, 3 + n] += BL11a_hexa[npoint][m, 3 + p] * l_perisp[p, n];
                            BL11[npoint][m, 6 + n] += BL11a_hexa[npoint][m, 6 + p] * l_perisp[p, n];
                        }
                    }
                }

                //
                for (int m = 0; m < 6; m++)
                {
                    for (int n = 0; n < 9; n++)
                    {
                        //BL1112sun01_hexa[npoint][m, n] = 0; sp1
                        for (int p = 0; p < 9; p++)
                        {
                            BL1112sun01_hexa[npoint][m, n] += BL11[npoint][m, p] * BL12_hexa[npoint][p, n];
                        }
                    }
                }
                for (int m = 0; m < 6; m++)
                {
                    for (int n = 0; n < 9; n++)
                    {
                        BL1112sun01_hexa[npoint][m, n] += BL01_hexa[npoint][m, n];
                    }
                }

                //
                for (int m = 0; m < 6; m++)
                {
                    for (int n = 0; n < 24; n++)
                    {
                        BL[npoint][m, n] = 0;
                        for (int p = 0; p < 9; p++)
                        {
                            BL[npoint][m, n] += BL1112sun01_hexa[npoint][m, p] * BL13_hexa[npoint][p, n];
                        }
                    }
                }

                //
                for (int m = 0; m < 24; m++)
                {
                    //fxk1[npoint][m] = 0; //sp1
                    for (int n = 0; n < 6; n++)
                    {
                        fxk1[npoint][m] += BL[npoint][n, m] * sunt_ol_Spkvec[npoint][n];
                    }
                }
            }

            //
            //for (int m = 0; m < 24; m++)
            //{
            //    fxk1[nGaussPoints][m] = 0; //sp1
            //}
            for (int npoint = 0; npoint < nGaussPoints; npoint++)
            {
                for (int m = 0; m < 24; m++)
                {
                    fxk1[nGaussPoints][m] += fxk1[npoint][m];
                }
            }

            return fxk1[nGaussPoints];
        }

        private void InitializeBland_sunt_ol_Spkvec() //.first_calc_for_Kmatrices
        {
            // DEN XREIAZETAI PIA INITIALIZE AFTWN pou kratiountai ginontai initialize ekei pou upologizontai 

            //double[][,] BL11 = new double[nGaussPoints][,];
            //double[][,] BL1112sun01_hexa = new double[nGaussPoints][,];
            //for (int npoint = 0; npoint < nGaussPoints; npoint++)
            //{
            //    BL11[npoint] = new double[6, 9];
            //    BL1112sun01_hexa[npoint] = new double[6, 9];
            //}

            for (int npoint = 0; npoint < nGaussPoints; npoint++)
            {

                //
                //for (int m = 0; m < 6; m++)
                //{
                //    sunt_ol_Spkvec[npoint][m] = sunt_oloklhrwmatos[npoint] * Spkvec[npoint][m];
                //}

                //
                //double[,] l_perisp = new double[3, 3];
                //for (int m = 0; m < 3; m++)
                //{
                //    for (int n = 0; n < 3; n++)
                //    {
                //        //l_perisp[m, n] = 0; // sp1
                //        for (int p = 0; p < 8; p++)
                //        {
                //            l_perisp[m, n] += ll1_hexa[npoint][m, p] * ll2[p, n];
                //        }
                //    }
                //}

                //
                //for (int m = 0; m < 6; m++)
                //{
                //    for (int n = 0; n < 9; n++)
                //    {
                //        BL11[npoint][m, n] = 0;
                //    }
                //} //sp1
                //for (int m = 0; m < 6; m++)
                //{
                //    for (int n = 0; n < 3; n++)
                //    {
                //        for (int p = 0; p < 3; p++)
                //        {
                //            BL11[npoint][m, n] += BL11a_hexa[npoint][m, p] * l_perisp[p, n];
                //            BL11[npoint][m, 3 + n] += BL11a_hexa[npoint][m, 3 + p] * l_perisp[p, n];
                //            BL11[npoint][m, 6 + n] += BL11a_hexa[npoint][m, 6 + p] * l_perisp[p, n];
                //        }
                //    }
                //}

                ////
                //for (int m = 0; m < 6; m++)
                //{
                //    for (int n = 0; n < 9; n++)
                //    {
                //        //BL1112sun01_hexa[npoint][m, n] = 0; //sp1
                //        for (int p = 0; p < 9; p++)
                //        {
                //            BL1112sun01_hexa[npoint][m, n] += BL11[npoint][m, p] * BL12_hexa[npoint][p, n];
                //        }
                //    }
                //}
                //for (int m = 0; m < 6; m++)
                //{
                //    for (int n = 0; n < 9; n++)
                //    {
                //        BL1112sun01_hexa[npoint][m, n] += BL01_hexa[npoint][m, n];
                //    }
                //}

                //
                //for (int m = 0; m < 6; m++)
                //{
                //    for (int n = 0; n < 24; n++)
                //    {
                //        BL[npoint][m, n] = 0;
                //        for (int p = 0; p < 9; p++)
                //        {
                //            BL[npoint][m, n] += BL1112sun01_hexa[npoint][m, p] * BL13_hexa[npoint][p, n];
                //        }
                //    }
                //}
            }

        }


        private double[,] UpdateKmatrices(IElement element)
        {
            double[,] k_stoixeiou = new double[24, 24];


            // KAI INITIALIZE TWN APARAITHTWN APO TIS FORCES pou den tha pothikevontai poia
            double[][] sunt_ol_Spkvec = new double[nGaussPoints][];
            double[][,] BL = new double[nGaussPoints][,];
            for (int gpoint = 0; gpoint < nGaussPoints; gpoint++)
            {
                sunt_ol_Spkvec[gpoint] = new double[6];
                BL[gpoint] = new double[6, 24];

            }

            //PRWTA TA APARAITHTA APO TIS FORCES POU DEN THA KRATIOUNTAI pia apo tis FORCES
            // upologismos entos forces olwn twn apaitoumenwn mhtrwwn
            double[,] ll2 = new double[8, 3];
            for (int m = 0; m < 8; m++)
            {
                for (int n = 0; n < 3; n++)
                {
                    ll2[m, n] = tu_i[m][n];
                }
            }
            // upologismos ennoeitai kai twn mhtrwwn apo tous arxikous upologismous pou xreiazontai edw (ola ews ll1_hexa)
            double[,] Ni_ksi;
            double[,] Ni_heta;
            double[,] Ni_zeta;
            Tuple<double[,], double[,], double[,], double[]> shapeFunctionDerivativesData;
            shapeFunctionDerivativesData = GetShapeFunctionDerivativesData(gp_d1_disp, gp_d2_disp, gp_d3_disp);
            Ni_ksi = shapeFunctionDerivativesData.Item1;
            Ni_heta = shapeFunctionDerivativesData.Item2;
            Ni_zeta = shapeFunctionDerivativesData.Item3;
            double[][,] ll1_hexa;
            ll1_hexa = Getll1Hexa(Ni_ksi, Ni_heta, Ni_zeta, gp_d1_disp, gp_d2_disp, gp_d3_disp);
            //11a 12 13 kai 01 epishs xreiazontai opote kai to J_0inv_hexa pou afta theloun 
            Tuple<double[][,], double[]> J_0inv_hexaAndDetJ_0;
            J_0inv_hexaAndDetJ_0 = GetJ_0invHexaAndDetJ_0(ll1_hexa, element);
            double[][,] J_0inv_hexa;
            J_0inv_hexa = J_0inv_hexaAndDetJ_0.Item1;
            double[][,] BL13_hexa;
            BL13_hexa = GetBL13Hexa(ll1_hexa, gp_d1_disp, gp_d2_disp, gp_d3_disp);
            double[][,] BL11a_hexa; // exoume tosa [,] osa einai kai ta gpoints
            double[][,] BL12_hexa;
            double[][,] BL01_hexa;
            BL11a_hexa = GetBL11a_hexa(J_0inv_hexa);
            BL12_hexa = GetBL12_hexa(J_0inv_hexa);
            BL01_hexa = GetBL01_hexa(J_0inv_hexa);

            double[][,] BL11 = new double[nGaussPoints][,];
            double[][,] BL1112sun01_hexa = new double[nGaussPoints][,];
            for (int npoint = 0; npoint < nGaussPoints; npoint++)
            {
                BL11[npoint] = new double[6, 9];
                BL1112sun01_hexa[npoint] = new double[6, 9];
            }



            for (int npoint = 0; npoint < nGaussPoints; npoint++)
            {

                //
                for (int m = 0; m < 6; m++)
                {
                    sunt_ol_Spkvec[npoint][m] = sunt_oloklhrwmatos[npoint] * Spkvec[npoint][m];
                }

                //
                double[,] l_perisp = new double[3, 3];
                for (int m = 0; m < 3; m++)
                {
                    for (int n = 0; n < 3; n++)
                    {
                        //l_perisp[m, n] = 0; //sp1
                        for (int p = 0; p < 8; p++)
                        {
                            l_perisp[m, n] += ll1_hexa[npoint][m, p] * ll2[p, n];
                        }
                    }
                }

                //
                //for (int m = 0; m < 6; m++)
                //{
                //    for (int n = 0; n < 9; n++)
                //    {
                //        BL11[npoint][m,n] = 0;
                //    }
                //} sp1
                for (int m = 0; m < 6; m++)
                {
                    for (int n = 0; n < 3; n++)
                    {
                        for (int p = 0; p < 3; p++)
                        {
                            BL11[npoint][m, n] += BL11a_hexa[npoint][m, p] * l_perisp[p, n];
                            BL11[npoint][m, 3 + n] += BL11a_hexa[npoint][m, 3 + p] * l_perisp[p, n];
                            BL11[npoint][m, 6 + n] += BL11a_hexa[npoint][m, 6 + p] * l_perisp[p, n];
                        }
                    }
                }

                //
                for (int m = 0; m < 6; m++)
                {
                    for (int n = 0; n < 9; n++)
                    {
                        //BL1112sun01_hexa[npoint][m, n] = 0; sp1
                        for (int p = 0; p < 9; p++)
                        {
                            BL1112sun01_hexa[npoint][m, n] += BL11[npoint][m, p] * BL12_hexa[npoint][p, n];
                        }
                    }
                }
                for (int m = 0; m < 6; m++)
                {
                    for (int n = 0; n < 9; n++)
                    {
                        BL1112sun01_hexa[npoint][m, n] += BL01_hexa[npoint][m, n];
                    }
                }

                //
                for (int m = 0; m < 6; m++)
                {
                    for (int n = 0; n < 24; n++)
                    {
                        BL[npoint][m, n] = 0;
                        for (int p = 0; p < 9; p++)
                        {
                            BL[npoint][m, n] += BL1112sun01_hexa[npoint][m, p] * BL13_hexa[npoint][p, n];
                        }
                    }
                }
            }
            //PRWTA TA APARAITHTA APO TIS FORCES POU DEN THA KRATIOUNTAI pia apo tis FORCES

            //DEFTERON UPOLOGIZETAI TO BNL pou den tha proupologizetai pleon apo to initial configuration
            double[][,] BNL1_hexa;
            double[][,] BNL_hexa;
            BNL1_hexa = GetBNL1_hexa(J_0inv_hexa);
            BNL_hexa = new double[nGaussPoints][,];
            for (int gpoint = 0; gpoint < nGaussPoints; gpoint++)
            {
                // initialize diastaseis twn mhtrwwn kai meta gemisma keliwn (olwn h mono oswn mporoume sthn arxh)
                BNL_hexa[gpoint] = new double[9, 24];

                //
                for (int m = 0; m < 9; m++)
                {
                    for (int n = 0; n < 24; n++)
                    {
                        //BNL_hexa[gpoint][m, n] = 0;    //sp1
                        for (int p = 0; p < 9; p++)
                        {
                            BNL_hexa[gpoint][m, n] += BNL1_hexa[gpoint][m, p] * BL13_hexa[gpoint][p, n];
                        }
                    }
                }
            }
            //DEFTERON UPOLOGIZETAI TO BNL pou den tha proupologizetai pleon apo to initial configuration





            double[,] sunt_ol_SPK_epi_BNL_hexa = new double[9, 24];
            double[,] sunt_ol_cons_disp = new double[6, 6];
            double[,] sunt_ol_cons_disp_epi_BL = new double[6, 24];

            double[][,] sunt_ol_Spk = new double[nGaussPoints][,];
            for (int npoint = 0; npoint < nGaussPoints; npoint++)
            {
                sunt_ol_Spk[npoint] = new double[3, 3];
            }

            double[][,] kl_ = new double[nGaussPoints + 1][,];
            double[][,] knl_ = new double[nGaussPoints + 1][,];
            for (int npoint = 0; npoint < nGaussPoints + 1; npoint++)
            {
                kl_[npoint] = new double[24, 24];
                knl_[npoint] = new double[24, 24];
            }



            for (int npoint = 0; npoint < nGaussPoints; npoint++)
            {
                //
                sunt_ol_Spk[npoint][0, 0] = sunt_ol_Spkvec[npoint][0];
                sunt_ol_Spk[npoint][0, 1] = sunt_ol_Spkvec[npoint][3];
                sunt_ol_Spk[npoint][0, 2] = sunt_ol_Spkvec[npoint][5];
                sunt_ol_Spk[npoint][1, 0] = sunt_ol_Spkvec[npoint][3];
                sunt_ol_Spk[npoint][1, 1] = sunt_ol_Spkvec[npoint][1];
                sunt_ol_Spk[npoint][1, 2] = sunt_ol_Spkvec[npoint][4];
                sunt_ol_Spk[npoint][2, 0] = sunt_ol_Spkvec[npoint][5];
                sunt_ol_Spk[npoint][2, 1] = sunt_ol_Spkvec[npoint][4];
                sunt_ol_Spk[npoint][2, 2] = sunt_ol_Spkvec[npoint][2];

                //
                for (int m = 0; m < 6; m++)
                {
                    for (int n = 0; n < 6; n++)
                    {
                        sunt_ol_cons_disp[m, n] = sunt_oloklhrwmatos[npoint] * Cons_disp[npoint][m, n];
                    }
                }

                //
                for (int m = 0; m < 6; m++)
                {
                    for (int n = 0; n < 24; n++)
                    {
                        sunt_ol_cons_disp_epi_BL[m, n] = 0;
                        for (int p = 0; p < 6; p++)
                        {
                            sunt_ol_cons_disp_epi_BL[m, n] += sunt_ol_cons_disp[m, p] * BL[npoint][p, n];
                        }
                    }
                }

                //
                for (int m = 0; m < 24; m++)
                {
                    for (int n = 0; n < 24; n++)
                    {
                        kl_[npoint][m, n] = 0;
                        for (int p = 0; p < 6; p++)
                        {
                            kl_[npoint][m, n] += BL[npoint][p, m] * sunt_ol_cons_disp_epi_BL[p, n];
                        }
                    }
                }
                //tha athroisoume meta ola ta kl- sthn teleftaia thesi

                //
                for (int m = 0; m < 3; m++) //prwtes 3x24 grammes
                {
                    for (int n = 0; n < 24; n++)
                    {
                        sunt_ol_SPK_epi_BNL_hexa[m, n] = 0;
                        sunt_ol_SPK_epi_BNL_hexa[3 + m, n] = 0;
                        sunt_ol_SPK_epi_BNL_hexa[6 + m, n] = 0;
                        for (int p = 0; p < 3; p++)
                        {
                            sunt_ol_SPK_epi_BNL_hexa[m, n] += sunt_ol_Spk[npoint][m, p] * BNL_hexa[npoint][p, n];
                            sunt_ol_SPK_epi_BNL_hexa[3 + m, n] += sunt_ol_Spk[npoint][m, p] * BNL_hexa[npoint][3 + p, n];
                            sunt_ol_SPK_epi_BNL_hexa[6 + m, n] += sunt_ol_Spk[npoint][m, p] * BNL_hexa[npoint][6 + p, n];
                        }
                    }
                }

                //
                for (int m = 0; m < 24; m++)
                {
                    for (int n = 0; n < 24; n++)
                    {
                        knl_[npoint][m, n] = 0;
                        for (int p = 0; p < 9; p++)
                        {
                            knl_[npoint][m, n] += BNL_hexa[npoint][p, m] * sunt_ol_SPK_epi_BNL_hexa[p, n];
                        }
                    }
                }
                //tha athroisoume meta ola ta knl_ sthn teleftaia thesi i kateftheian sto k_stoixeiou

            }

            // athroisma olwn twn gpoints se k_stoixeiou kai prwta mhdenismos aftou
            for (int m = 0; m < 24; m++)
            {
                for (int n = 0; n < 24; n++)
                {
                    kl_[nGaussPoints][m, n] = 0;
                    knl_[nGaussPoints][m, n] = 0;
                }
            }
            for (int npoint = 0; npoint < nGaussPoints; npoint++)
            {
                for (int m = 0; m < 24; m++)
                {
                    for (int n = 0; n < 24; n++)
                    {
                        kl_[nGaussPoints][m, n] += kl_[npoint][m, n];
                        knl_[nGaussPoints][m, n] += knl_[npoint][m, n];
                    }
                }
            }
            for (int m = 0; m < 24; m++)
            {
                for (int n = 0; n < 24; n++)
                {
                    k_stoixeiou[m, n] = kl_[nGaussPoints][m, n] + knl_[nGaussPoints][m, n];
                }
            }

            return k_stoixeiou;
        }

        // telikes entoles kai mhtrwo mazas apo to hexa8

        public Tuple<double[], double[]> CalculateStresses(Element element, double[] localTotalDisplacements, double[] localdDisplacements)
        {
            this.UpdateCoordinateData(localTotalDisplacements, element);
            for (int npoint = 0; npoint < materialsAtGaussPoints.Length; npoint++)
            {
                materialsAtGaussPoints[npoint].UpdateMaterial(new StressStrainVectorContinuum3D(GLvec[npoint]));
            }
            return new Tuple<double[], double[]>(GLvec[materialsAtGaussPoints.Length - 1], materialsAtGaussPoints[materialsAtGaussPoints.Length - 1].Stresses.Data);
            //TODO mono to teleftaio dianusma tha epistrefei?
        }

        public double[] CalculateForces(Element element, double[] localTotalDisplacements, double[] localdDisplacements)
        {
            double[] fxk1;
            for (int npoint = 0; npoint < materialsAtGaussPoints.Length; npoint++)
            {
                for (int j = 0; j < 6; j++)
                { Spkvec[npoint][j] = materialsAtGaussPoints[npoint].Stresses[j]; }
            }
            fxk1 = this.UpdateForces(element);
            return fxk1;
        }

        public double[] CalculateForcesForLogging(Element element, double[] localDisplacements)
        {
            return CalculateForces(element, localDisplacements, new double[localDisplacements.Length]);
        }

        public virtual IMatrix2D StiffnessMatrix(IElement element)
        {
            double[,] k_stoixeiou = new double[24, 24];
            if (Cons_disp == null)
            {
                //this.CalculateShapeFunctionAndGaussPointData();
                //this.GetInitialGeometricDataAndInitializeMatrices(element);
                this.CalculateInitialConfigurationData(element);

                this.UpdateCoordinateData(new double[24], element);
                for (int npoint = 0; npoint < materialsAtGaussPoints.Length; npoint++)// loop gia getfirstStressesFromMaterial
                {
                    for (int j = 0; j < 6; j++)
                    { Spkvec[npoint][j] = materialsAtGaussPoints[npoint].Stresses[j]; }
                }
                this.InitializeBland_sunt_ol_Spkvec();// meta to get twn stresses apo to material dioiti periexei ton pol/smo suntol epi Spkvec
            }
            for (int npoint = 0; npoint < materialsAtGaussPoints.Length; npoint++)
            {
                for (int j = 0; j < 6; j++)
                {
                    for (int k = 0; k < 6; k++)
                    { Cons_disp[npoint][j, k] = materialsAtGaussPoints[npoint].ConstitutiveMatrix[j, k]; }
                }
            }
            k_stoixeiou = this.UpdateKmatrices(element);
            IMatrix2D element_stiffnessMatrix = new Matrix2D(k_stoixeiou); // TODO giati de ginetai return dof.Enumerator.GetTransformedMatrix, xrhsh symmetric
            return element_stiffnessMatrix;
        }

        public bool MaterialModified
        {
            get
            {
                foreach (IContinuumMaterial3D material in materialsAtGaussPoints)
                    if (material.Modified) return true;
                return false;
            }
        }

        public void ResetMaterialModified()
        {
            foreach (IContinuumMaterial3D material in materialsAtGaussPoints) material.ResetModified();
        }

        public void ClearMaterialState()
        {
            foreach (IContinuumMaterial3D m in materialsAtGaussPoints) m.ClearState();
        }

        public void SaveMaterialState()
        {
            foreach (IContinuumMaterial3D m in materialsAtGaussPoints) m.SaveState();
        }

        public void ClearMaterialStresses()
        {
            foreach (IContinuumMaterial3D m in materialsAtGaussPoints) m.ClearStresses();
        }

        // omoiws me hexa 8 shell8disp implemented
        public int ID
        {
            get { return 13; }
        }
        public ElementDimensions ElementDimensions
        {
            get { return ElementDimensions.ThreeD; }
        }

        public IElementDOFEnumerator DOFEnumerator
        {
            get { return dofEnumerator; }
            set { dofEnumerator = value; }
        }

        public virtual IList<IList<DOFType>> GetElementDOFTypes(IElement element)
        {
            return dofTypes;
        }

        //NOT IMPLEMENTED
        public double[] CalculateAccelerationForces(Element element, IList<MassAccelerationLoad> loads)
        {
            return new double[24];
        }

        public virtual IMatrix2D MassMatrix(IElement element)
        {
            return new Matrix2D(24, 24);
        }

        public virtual IMatrix2D DampingMatrix(IElement element)
        {

            return new Matrix2D(24, 24);
        }
        //NOT IMPLEMENTED ews edw


        // Perioxh IEmbeddedHostElement 

        // methodoi apo to Hexa8.cs
        protected double[,] GetCoordinatesTranspose(Element element)
        {
            double[,] faXYZ = new double[3, dofTypes.Length];
            for (int i = 0; i < dofTypes.Length; i++)
            {
                faXYZ[0, i] = element.Nodes[i].X;
                faXYZ[1, i] = element.Nodes[i].Y;
                faXYZ[2, i] = element.Nodes[i].Z;
            }
            return faXYZ;
        }

        private double[] CalcH8Shape(double fXi, double fEta, double fZeta)
        {
            const double fSqC125 = 0.5;
            double fXiP = (1.0 + fXi) * fSqC125;
            double fEtaP = (1.0 + fEta) * fSqC125;
            double fZetaP = (1.0 + fZeta) * fSqC125;
            double fXiM = (1.0 - fXi) * fSqC125;
            double fEtaM = (1.0 - fEta) * fSqC125;
            double fZetaM = (1.0 - fZeta) * fSqC125;

            double[] H8Shape = new double[8]; // PROSTHIKI EMBEDDED : allagh twn Calc8shape pou eginan copy

            H8Shape[6] = fXiM * fEtaM * fZetaM;
            H8Shape[7] = fXiP * fEtaM * fZetaM;
            H8Shape[4] = fXiP * fEtaP * fZetaM;
            H8Shape[5] = fXiM * fEtaP * fZetaM;
            H8Shape[2] = fXiM * fEtaM * fZetaP;
            H8Shape[3] = fXiP * fEtaM * fZetaP;
            H8Shape[0] = fXiP * fEtaP * fZetaP;
            H8Shape[1] = fXiM * fEtaP * fZetaP;

            return H8Shape;

            //return new double[]
            //{
            //    fXiM * fEtaM * fZetaM,
            //    fXiP * fEtaM * fZetaM,
            //    fXiP * fEtaP * fZetaM,
            //    fXiM * fEtaP * fZetaM,
            //    fXiM * fEtaM * fZetaP,
            //    fXiP * fEtaM * fZetaP,
            //    fXiP * fEtaP * fZetaP,
            //    fXiM * fEtaP * fZetaP
            //};
        }

        private double[] CalcH8NablaShape(double fXi, double fEta, double fZeta)
        {
            const double fSq125 = 0.35355339059327376220042218105242;
            double fXiP = (1.0 + fXi) * fSq125;
            double fEtaP = (1.0 + fEta) * fSq125;
            double fZetaP = (1.0 + fZeta) * fSq125;
            double fXiM = (1.0 - fXi) * fSq125;
            double fEtaM = (1.0 - fEta) * fSq125;
            double fZetaM = (1.0 - fZeta) * fSq125;

            double[] faDS = new double[24];
            faDS[6] = -fEtaM * fZetaM;
            faDS[4] = fEtaP * fZetaM;
            faDS[2] = -fEtaM * fZetaP;
            faDS[0] = fEtaP * fZetaP;
            faDS[7] = -faDS[6];
            faDS[5] = -faDS[4];
            faDS[3] = -faDS[2];
            faDS[1] = -faDS[0];


            faDS[14] = -fXiM * fZetaM;
            faDS[15] = -fXiP * fZetaM;
            faDS[10] = -fXiM * fZetaP;
            faDS[11] = -fXiP * fZetaP;
            faDS[12] = -faDS[15];
            faDS[13] = -faDS[14];
            faDS[8] = -faDS[11];
            faDS[9] = -faDS[10];


            faDS[22] = -fXiM * fEtaM;
            faDS[23] = -fXiP * fEtaM;
            faDS[20] = -fXiP * fEtaP;
            faDS[21] = -fXiM * fEtaP;
            faDS[18] = -faDS[22];
            faDS[19] = -faDS[23];
            faDS[16] = -faDS[20];
            faDS[17] = -faDS[21];

            return faDS;
        }

        protected static double determinantTolerance = 0.00000001;
        private Tuple<double[,], double[,], double> CalcH8JDetJ(double[,] faXYZ, double[] faDS)
        {
            double[,] faJ = new double[3, 3];
            faJ[0, 0] = faDS[0] * faXYZ[0, 0] + faDS[1] * faXYZ[0, 1] + faDS[2] * faXYZ[0, 2] + faDS[3] * faXYZ[0, 3] + faDS[4] * faXYZ[0, 4] + faDS[5] * faXYZ[0, 5] + faDS[6] * faXYZ[0, 6] + faDS[7] * faXYZ[0, 7];
            faJ[0, 1] = faDS[0] * faXYZ[1, 0] + faDS[1] * faXYZ[1, 1] + faDS[2] * faXYZ[1, 2] + faDS[3] * faXYZ[1, 3] + faDS[4] * faXYZ[1, 4] + faDS[5] * faXYZ[1, 5] + faDS[6] * faXYZ[1, 6] + faDS[7] * faXYZ[1, 7];
            faJ[0, 2] = faDS[0] * faXYZ[2, 0] + faDS[1] * faXYZ[2, 1] + faDS[2] * faXYZ[2, 2] + faDS[3] * faXYZ[2, 3] + faDS[4] * faXYZ[2, 4] + faDS[5] * faXYZ[2, 5] + faDS[6] * faXYZ[2, 6] + faDS[7] * faXYZ[2, 7];
            faJ[1, 0] = faDS[8] * faXYZ[0, 0] + faDS[9] * faXYZ[0, 1] + faDS[10] * faXYZ[0, 2] + faDS[11] * faXYZ[0, 3] + faDS[12] * faXYZ[0, 4] + faDS[13] * faXYZ[0, 5] + faDS[14] * faXYZ[0, 6] + faDS[15] * faXYZ[0, 7];
            faJ[1, 1] = faDS[8] * faXYZ[1, 0] + faDS[9] * faXYZ[1, 1] + faDS[10] * faXYZ[1, 2] + faDS[11] * faXYZ[1, 3] + faDS[12] * faXYZ[1, 4] + faDS[13] * faXYZ[1, 5] + faDS[14] * faXYZ[1, 6] + faDS[15] * faXYZ[1, 7];
            faJ[1, 2] = faDS[8] * faXYZ[2, 0] + faDS[9] * faXYZ[2, 1] + faDS[10] * faXYZ[2, 2] + faDS[11] * faXYZ[2, 3] + faDS[12] * faXYZ[2, 4] + faDS[13] * faXYZ[2, 5] + faDS[14] * faXYZ[2, 6] + faDS[15] * faXYZ[2, 7];
            faJ[2, 0] = faDS[16] * faXYZ[0, 0] + faDS[17] * faXYZ[0, 1] + faDS[18] * faXYZ[0, 2] + faDS[19] * faXYZ[0, 3] + faDS[20] * faXYZ[0, 4] + faDS[21] * faXYZ[0, 5] + faDS[22] * faXYZ[0, 6] + faDS[23] * faXYZ[0, 7];
            faJ[2, 1] = faDS[16] * faXYZ[1, 0] + faDS[17] * faXYZ[1, 1] + faDS[18] * faXYZ[1, 2] + faDS[19] * faXYZ[1, 3] + faDS[20] * faXYZ[1, 4] + faDS[21] * faXYZ[1, 5] + faDS[22] * faXYZ[1, 6] + faDS[23] * faXYZ[1, 7];
            faJ[2, 2] = faDS[16] * faXYZ[2, 0] + faDS[17] * faXYZ[2, 1] + faDS[18] * faXYZ[2, 2] + faDS[19] * faXYZ[2, 3] + faDS[20] * faXYZ[2, 4] + faDS[21] * faXYZ[2, 5] + faDS[22] * faXYZ[2, 6] + faDS[23] * faXYZ[2, 7];

            double fDet1 = faJ[0, 0] * (faJ[1, 1] * faJ[2, 2] - faJ[2, 1] * faJ[1, 2]);
            double fDet2 = -faJ[0, 1] * (faJ[1, 0] * faJ[2, 2] - faJ[2, 0] * faJ[1, 2]);
            double fDet3 = faJ[0, 2] * (faJ[1, 0] * faJ[2, 1] - faJ[2, 0] * faJ[1, 1]);
            double fDetJ = fDet1 + fDet2 + fDet3;
            if (fDetJ < determinantTolerance)
                throw new ArgumentException(String.Format("Jacobian determinant is negative or under tolerance ({0} < {1}). Check the order of nodes or the element geometry.", fDetJ, determinantTolerance));

            double fDetInv = 1.0 / fDetJ;
            double[,] faJInv = new double[3, 3];
            faJInv[0, 0] = (faJ[1, 1] * faJ[2, 2] - faJ[2, 1] * faJ[1, 2]) * fDetInv;
            faJInv[1, 0] = (faJ[2, 0] * faJ[1, 2] - faJ[1, 0] * faJ[2, 2]) * fDetInv;
            faJInv[2, 0] = (faJ[1, 0] * faJ[2, 1] - faJ[2, 0] * faJ[1, 1]) * fDetInv;
            faJInv[0, 1] = (faJ[2, 1] * faJ[0, 2] - faJ[0, 1] * faJ[2, 2]) * fDetInv;
            faJInv[1, 1] = (faJ[0, 0] * faJ[2, 2] - faJ[2, 0] * faJ[0, 2]) * fDetInv;
            faJInv[2, 1] = (faJ[2, 0] * faJ[0, 1] - faJ[2, 1] * faJ[0, 0]) * fDetInv;
            faJInv[0, 2] = (faJ[0, 1] * faJ[1, 2] - faJ[1, 1] * faJ[0, 2]) * fDetInv;
            faJInv[1, 2] = (faJ[1, 0] * faJ[0, 2] - faJ[0, 0] * faJ[1, 2]) * fDetInv;
            faJInv[2, 2] = (faJ[0, 0] * faJ[1, 1] - faJ[1, 0] * faJ[0, 1]) * fDetInv;

            return new Tuple<double[,], double[,], double>(faJ, faJInv, fDetJ);
        }

        public EmbeddedNode BuildHostElementEmbeddedNode(Element element, Node node, IEmbeddedDOFInHostTransformationVector transformationVector)
        {
            var points = GetNaturalCoordinates(element, node);
            if (points.Length == 0) return null;

            element.EmbeddedNodes.Add(node);
            var embeddedNode = new EmbeddedNode(node, element, transformationVector.GetDependentDOFTypes);
            for (int i = 0; i < points.Length; i++)
                embeddedNode.Coordinates.Add(points[i]);
            return embeddedNode;
        }

        private double[] GetNaturalCoordinates(Element element, Node node)
        {
            double[] mins = new double[] { element.Nodes[0].X, element.Nodes[0].Y, element.Nodes[0].Z };
            double[] maxes = new double[] { element.Nodes[0].X, element.Nodes[0].Y, element.Nodes[0].Z };
            for (int i = 0; i < element.Nodes.Count; i++)
            {
                mins[0] = mins[0] > element.Nodes[i].X ? element.Nodes[i].X : mins[0];
                mins[1] = mins[1] > element.Nodes[i].Y ? element.Nodes[i].Y : mins[1];
                mins[2] = mins[2] > element.Nodes[i].Z ? element.Nodes[i].Z : mins[2];
                maxes[0] = maxes[0] < element.Nodes[i].X ? element.Nodes[i].X : maxes[0];
                maxes[1] = maxes[1] < element.Nodes[i].Y ? element.Nodes[i].Y : maxes[1];
                maxes[2] = maxes[2] < element.Nodes[i].Z ? element.Nodes[i].Z : maxes[2];
            }
            //return new double[] { (node.X - mins[0]) / ((maxes[0] - mins[0]) / 2) - 1,
            //    (node.Y - mins[1]) / ((maxes[1] - mins[1]) / 2) - 1,
            //    (node.Z - mins[2]) / ((maxes[2] - mins[2]) / 2) - 1 };

            bool maybeInsideElement = node.X <= maxes[0] && node.X >= mins[0] &&
                node.Y <= maxes[1] && node.Y >= mins[1] &&
                node.Z <= maxes[2] && node.Z >= mins[2];
            if (maybeInsideElement == false) return new double[0];

            const int jacobianSize = 3;
            const int maxIterations = 1000;
            const double tolerance = 1e-10;
            int iterations = 0;
            double deltaNaturalCoordinatesNormSquare = 100;
            double[] naturalCoordinates = new double[] { 0, 0, 0 };
            const double toleranceSquare = tolerance * tolerance;

            while (deltaNaturalCoordinatesNormSquare > toleranceSquare && iterations < maxIterations)
            {
                iterations++;
                var shapeFunctions = CalcH8Shape(naturalCoordinates[0], naturalCoordinates[1], naturalCoordinates[2]);
                double[] coordinateDifferences = new double[] { 0, 0, 0 };
                for (int i = 0; i < shapeFunctions.Length; i++)
                {
                    coordinateDifferences[0] += shapeFunctions[i] * element.Nodes[i].X;
                    coordinateDifferences[1] += shapeFunctions[i] * element.Nodes[i].Y;
                    coordinateDifferences[2] += shapeFunctions[i] * element.Nodes[i].Z;
                }
                coordinateDifferences[0] = node.X - coordinateDifferences[0];
                coordinateDifferences[1] = node.Y - coordinateDifferences[1];
                coordinateDifferences[2] = node.Z - coordinateDifferences[2];

                double[,] faXYZ = GetCoordinatesTranspose(element);
                double[] nablaShapeFunctions = CalcH8NablaShape(naturalCoordinates[0], naturalCoordinates[1], naturalCoordinates[2]);
                var inverseJacobian = CalcH8JDetJ(faXYZ, nablaShapeFunctions).Item2;

                double[] deltaNaturalCoordinates = new double[] { 0, 0, 0 };
                for (int i = 0; i < jacobianSize; i++)
                    for (int j = 0; j < jacobianSize; j++)
                        deltaNaturalCoordinates[i] += inverseJacobian[j, i] * coordinateDifferences[j];
                for (int i = 0; i < 3; i++)
                    naturalCoordinates[i] += deltaNaturalCoordinates[i];

                deltaNaturalCoordinatesNormSquare = 0;
                for (int i = 0; i < 3; i++)
                    deltaNaturalCoordinatesNormSquare += deltaNaturalCoordinates[i] * deltaNaturalCoordinates[i];
                //deltaNaturalCoordinatesNormSquare = Math.Sqrt(deltaNaturalCoordinatesNormSquare);
            }

            return naturalCoordinates.Count(x => Math.Abs(x) - 1.0 > tolerance) > 0 ? new double[0] : naturalCoordinates;
        }

        public double[] GetShapeFunctionsForNode(Element element, EmbeddedNode node)
        {
            double[,] elementCoordinates = GetCoordinatesTranspose(element);
            var shapeFunctions = CalcH8Shape(node.Coordinates[0], node.Coordinates[1], node.Coordinates[2]);
            var nablaShapeFunctions = CalcH8NablaShape(node.Coordinates[0], node.Coordinates[1], node.Coordinates[2]);
            var jacobian = CalcH8JDetJ(elementCoordinates, nablaShapeFunctions);

            return new double[]
            {
                shapeFunctions[0], shapeFunctions[1], shapeFunctions[2], shapeFunctions[3], shapeFunctions[4], shapeFunctions[5], shapeFunctions[6], shapeFunctions[7],
                nablaShapeFunctions[0], nablaShapeFunctions[1], nablaShapeFunctions[2], nablaShapeFunctions[3], nablaShapeFunctions[4], nablaShapeFunctions[5], nablaShapeFunctions[6], nablaShapeFunctions[7],
                nablaShapeFunctions[8], nablaShapeFunctions[9], nablaShapeFunctions[10], nablaShapeFunctions[11], nablaShapeFunctions[12], nablaShapeFunctions[13], nablaShapeFunctions[14], nablaShapeFunctions[15],
                nablaShapeFunctions[16], nablaShapeFunctions[17], nablaShapeFunctions[18], nablaShapeFunctions[19], nablaShapeFunctions[20], nablaShapeFunctions[21], nablaShapeFunctions[22], nablaShapeFunctions[23],
                jacobian.Item1[0, 0], jacobian.Item1[0, 1], jacobian.Item1[0, 2], jacobian.Item1[1, 0], jacobian.Item1[1, 1], jacobian.Item1[1, 2], jacobian.Item1[2, 0], jacobian.Item1[2, 1], jacobian.Item1[2, 2],
                jacobian.Item2[0, 0], jacobian.Item2[0, 1], jacobian.Item2[0, 2], jacobian.Item2[1, 0], jacobian.Item2[1, 1], jacobian.Item2[1, 2], jacobian.Item2[2, 0], jacobian.Item2[2, 1], jacobian.Item2[2, 2]
            };

            //double[] coords = GetNaturalCoordinates(element, node);
            //return CalcH8Shape(coords[0], coords[1], coords[2]);

            //double fXiP = (1.0 + coords[0]) * 0.5;
            //double fEtaP = (1.0 + coords[1]) * 0.5;
            //double fZetaP = (1.0 + coords[2]) * 0.5;
            //double fXiM = (1.0 - coords[0]) * 0.5;
            //double fEtaM = (1.0 - coords[1]) * 0.5;
            //double fZetaM = (1.0 - coords[2]) * 0.5;

            //return new double[] { fXiM * fEtaM * fZetaM,
            //    fXiP * fEtaM * fZetaM,
            //    fXiP * fEtaP * fZetaM,
            //    fXiM * fEtaP * fZetaM,
            //    fXiM * fEtaM * fZetaP,
            //    fXiP * fEtaM * fZetaP,
            //    fXiP * fEtaP * fZetaP,
            //    fXiM * fEtaP * fZetaP };
        }

        // Perioxh IEmbeddedHostElement

    }


}
