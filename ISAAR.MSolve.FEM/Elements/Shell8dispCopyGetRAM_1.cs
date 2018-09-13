using ISAAR.MSolve.Discretization;
using ISAAR.MSolve.Discretization.Interfaces;
//using ISAAR.MSolve.FEM.Embedding;//using ISAAR.MSolve.PreProcessor.Embedding;
// compa
using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.FEM.Interfaces;//using ISAAR.MSolve.PreProcessor.Interfaces;
using ISAAR.MSolve.FEM.Interpolation;
using ISAAR.MSolve.Materials.Interfaces;
using ISAAR.MSolve.Numerical.LinearAlgebra;//using ISAAR.MSolve.Matrices;
using ISAAR.MSolve.Numerical.LinearAlgebra.Interfaces;//using ISAAR.MSolve.Matrices.Interfaces;
using System;
using System.Collections.Generic;

namespace ISAAR.MSolve.FEM.Elements
{
    public class Shell8dispCopyGetRAM_1 : IStructuralFiniteElement
    {
        //metavlhtes opws sto hexa8
        protected readonly static DOFType[] nodalDOFTypes = new DOFType[] { DOFType.X, DOFType.Y, DOFType.Z, DOFType.RotX, DOFType.RotY };
        protected readonly static DOFType[][] dofTypes = new DOFType[][] { nodalDOFTypes, nodalDOFTypes, nodalDOFTypes,
            nodalDOFTypes, nodalDOFTypes, nodalDOFTypes, nodalDOFTypes, nodalDOFTypes };
        protected readonly IShellMaterial[] materialsAtGaussPoints; //compa : isotropic  //protected readonly IFiniteElementMaterial3D[] materialsAtGaussPoints;
        protected IElementDOFEnumerator dofEnumerator = new GenericDOFEnumerator();
        // ews edw 

        public double[][] oVn_i { get; set; }
        //public double[][] oV1_i { get; set; }
        //public double[][] oV2_i { get; set; }
        public int gp_d1 { get; set; } // den prepei na einai static--> shmainei idio gia ola taantikeimena afthw ths klashs
        public int gp_d2 { get; set; }
        public int gp_d3 { get; set; }
        public double[] tk { get; set; } //public static int[] tk { get; set; }
        private int nGaussPoints;

        //private double ksi;
        //private double heta;
        //private double zeta;
        private int npoint;

        private double[] integrationCoefficient;
        //private double[] a_123g;
        //private double a_1g;
        //private double a_2g;
        //private double a_3g;
        private readonly InterpolationShell8 interpolation;

        //private double[][] gausscoordinates;//3 dianysmata me tis timew tvn ksi heta zeta se ola ta gauss points
        //private double[][] shapeFunctions;// 8 dianusmata me tis times twn N1....N8 se kathe gauss point
        //private double[][] shapeFunctionDerivatives;// 16 dianusmata me tis times twn N1ksi....N8ksi,N1heta,....N8heta se kathe gauss point
        //private double[][,] ll1;//einai teliko kai oxi prok
        //private double[][,] J_0a;//einai teliko kai oxi prok
        // metavlhtes pou einai ex oloklhrou proupologismenes// PRWTA APO AFTO THA EKTELESTEI GETINTIALGEOMETRICDATA(ELEMENT)
        //private double[,] J_0b;    //einai idio gia ola ta gauss points 
        //private double[][,] J_0;       //den einai to idio gia ola ta gausspoint // einai teliko kai oxi prok
        //private double[] detJ_0; //[] osa kai ta gauss points
        //private double[][,] J_0inv;
        //private double[][,] BL11a;
        //private double[][,] BL12;
        //private double[][,] BNL1;

        //private double[][] tx_i; //8 arrays twn 3 stoixeiwn //den einai apo afta pou orizei o xrhsths
        private double[][] tU;   //8 arrays twn 6 stoixeiwn 
        private double[][] tUvec;//8 arrays twn 6 stoixeiwn

        //private double[][,] BL13; 
        //private double[] E;
        //private double[] ni;
        private double[][,] ConsCartes; // TODOGer delete

        private double[][] GLvec;
        //private double[][] SPKvec;



        public Shell8dispCopyGetRAM_1(IShellMaterial material, int gp_d1c, int gp_d2c, int gp_d3c) // compa isotropic
        {
            this.gp_d1 = gp_d1c;
            this.gp_d2 = gp_d2c;
            this.gp_d3 = gp_d3c;
            this.nGaussPoints = this.gp_d1 * this.gp_d2 * this.gp_d3;
            materialsAtGaussPoints = new IShellMaterial[nGaussPoints];
            for (int i = 0; i < nGaussPoints; i++) materialsAtGaussPoints[i] = material.Clone();
               
            this.interpolation = InterpolationShell8.UniqueInstance;
        }

        //public Shell8dispCopyGet(IFiniteElementMaterial3D material, IFiniteElementDOFEnumerator dofEnumerator)//pithanotata den xreiazetai
        //    : this(material, gp_d1, gp_d2, gp_d3)
        //{
        //    this.dofEnumerator = dofEnumerator;
        //}
        // ews edw


        

        private double[][,] GetBL11a(double [][,] J_0inv)
        {
            double[][,] BL11a;
            BL11a = new double[nGaussPoints][,];
            for (int j = 0; j < nGaussPoints; j++)
            { BL11a[j] = new double[6, 9]; }
            for (int j = 0; j < nGaussPoints; j++)
            {
                for (int k = 0; k < 6; k++)
                {
                    for (int l = 0; l < 9; l++)
                    { BL11a[j][k, l] = 0; }
                }

                for (int k = 0; k < 3; k++)
                {
                    for (int l = 0; l < 3; l++)
                    { BL11a[j][k, 3 * k + l] = J_0inv[j][k, l]; }
                }

                //gemisma [4,4] ews [4,6] kai [5,7] ews [5,9]
                for (int k = 0; k < 2; k++)
                {
                    for (int l = 0; l < 3; l++)
                    { BL11a[j][3 + k, 3 + 3 * k + l] = J_0inv[j][k, l]; }
                }

                //gemisma [4,1] ews [4,3] kai [5,4] ews [5,6]
                for (int k = 0; k < 2; k++)
                {
                    for (int l = 0; l < 3; l++)
                    { BL11a[j][3 + k, 3 * k + l] = J_0inv[j][1 + k, l]; }
                }

                for (int l = 0; l < 3; l++)
                { BL11a[j][5, l] = J_0inv[j][2, l]; }

                for (int l = 0; l < 3; l++)
                { BL11a[j][5, 6 + l] = J_0inv[j][0, l]; }
            }

            return BL11a;
        }

        private double[][,] GetBL12(double[][,] J_0inv)
        {
            double[][,] BL12;
            BL12 = new double[nGaussPoints][,];
            for (int j = 0; j < nGaussPoints; j++)
            { BL12[j] = new double[9, 9]; }
            for (int j = 0; j < nGaussPoints; j++)
            {
                for (int k = 0; k < 9; k++)
                {
                    for (int l = 0; l < 9; l++) BL12[j][k, l] = 0; 
                }

                for (int k = 0; k < 3; k++)
                {
                    for (int l = 0; l < 3; l++)
                    { BL12[j][k, 3 * k + l] = J_0inv[j][0, l]; }
                }

                for (int k = 0; k < 3; k++)
                {
                    for (int l = 0; l < 3; l++)
                    { BL12[j][3 + k, 3 * k + l] = J_0inv[j][1, l]; }
                }

                for (int k = 0; k < 3; k++)
                {
                    for (int l = 0; l < 3; l++)
                    { BL12[j][6 + k, 3 * k + l] = J_0inv[j][2, l]; }
                }
            }

            return BL12;
        }

        private double[][,] GetBL13(double[][] shapeFunctionDerivatives,double [][] tUvec , double[][,] J_0a)
        {
            double[][,] BL13;
            BL13 = new double[nGaussPoints][,];
            for (int j = 0; j < nGaussPoints; j++)
            { BL13[j] = new double[9, 40]; }
            for (int j = 0; j < nGaussPoints; j++)
            {
                //for (int k = 0; k < 9; k++)
                //{
                //    for (int l = 0; l < 40; l++)
                //    {
                //        BL13[j][k, l] = 0;
                //    }
                //} //sp1

                //sthles 1:3
                for (int m = 0; m < 8; m++)
                {
                    for (int k = 0; k < 3; k++)
                    {
                        for (int l = 0; l < 2; l++)
                        {
                            BL13[j][3 * k + l, 5 * m + k] = shapeFunctionDerivatives[m + 8 * l][j];
                        }
                    }
                }

                //sthles 4:5
                for (int m = 0; m < 8; m++)
                {
                    for (int k = 0; k < 3; k++)
                    {
                        for (int l = 0; l < 3; l++)
                        {
                            BL13[j][3 * k + l, 5 * m + 3] = -J_0a[j][l, m * 2 + 1] * tUvec[m][3 + k];
                            BL13[j][3 * k + l, 5 * m + 4] = +J_0a[j][l, m * 2 + 1] * tUvec[m][k];
                        }
                    }
                }
            }

            return BL13;
        }

        private double[][,] GetBNL1(double[][,] J_0inv)
        {
            double[][,] BNL1;
            BNL1 = new double[nGaussPoints][,];
            for (int j = 0; j < nGaussPoints; j++)
            { BNL1[j] = new double[9, 9]; }
            for (int j = 0; j < nGaussPoints; j++)
            {
                //for (int k = 0; k < 9; k++)
                //{
                //    for (int l = 0; l < 9; l++)
                //    { BNL1[j][k, l] = 0; }
                //} //sp1

                for (int m = 0; m < 3; m++)
                {
                    for (int k = 0; k < 3; k++)
                    {
                        for (int l = 0; l < 3; l++)
                        { BNL1[j][3 * m + k, 3 * m + l] = J_0inv[j][k, l]; }
                    }
                }
            }

            return BNL1;
        }

        private void CalculateInitialConfigurationData(IElement element, out double[][] tx_i)
        {
            // prosthiki ram 
            
            //double[,] J_0b;    //einai idio gia ola ta gauss points 
            //double[][,] J_0;       //den einai to idio gia ola ta gausspoint // einai teliko kai oxi prok
            double[] E;
            double[] ni;

            (double[][] gausscoordinates, double[][] shapeFunctions, double[][] shapeFunctionDerivatives, double[] a_123g) =
                interpolation.GetShapeFunctions(gp_d1, gp_d2, gp_d3);

            (double[][,] ll1, double[][,] J_0a) = JacobianShell8Calculations.Getll1AndJ_0a(
                nGaussPoints, tk, gausscoordinates, shapeFunctions, shapeFunctionDerivatives);
           
            double tV1norm;
            tx_i = new double[8][];
            tU = new double[8][];
            tUvec = new double[8][];
            double[][] oV1_i = new double[8][]; //tangent vector ''1'' initial configuration
            for (int j = 0; j < 8; j++)
            {
                tx_i[j] = new double[] { element.INodes[j].X, element.INodes[j].Y, element.INodes[j].Z, };
                tU[j] = new double[6];
                tUvec[j] = new double[6];
                oV1_i[j] = new double[3];
                for (int k = 0; k < 3; k++) { tU[j][3 + k] = oVn_i[j][k]; }

                tUvec[j][0] = tU[j][5];
                tUvec[j][1] = 0;
                tUvec[j][2] = -tU[j][3];

                tV1norm = Math.Sqrt(tUvec[j][0] * tUvec[j][0] + tUvec[j][1] * tUvec[j][1] + tUvec[j][2] * tUvec[j][2]);

                tUvec[j][0] = tUvec[j][0] / tV1norm;
                tUvec[j][1] = tUvec[j][1] / tV1norm;
                tUvec[j][2] = tUvec[j][2] / tV1norm;

                oV1_i[j][0] = tUvec[j][0];
                oV1_i[j][1] = tUvec[j][1];
                oV1_i[j][2] = tUvec[j][2];

                tUvec[j][3] = tU[j][3 + 1] * tUvec[j][2] - tU[j][3 + 2] * tUvec[j][1];
                tUvec[j][4] = tU[j][3 + 2] * tUvec[j][0] - tU[j][3 + 0] * tUvec[j][2];
                tUvec[j][5] = tU[j][3 + 0] * tUvec[j][1] - tU[j][3 + 1] * tUvec[j][0];
            }

            (double[][,] J_0inv, double[] detJ_0) =
                JacobianShell8Calculations.GetJ_0invAndDetJ_0(J_0a, element.INodes, oVn_i, nGaussPoints);


            //BL11a = GetBL11a(J_0inv);
            //BL12 = GetBL12(J_0inv);
            //BNL1 = GetBNL1(J_0inv);


            //PROSTHIKI gia ram 
            //double[,] Cons;
            //double[,] Cons_T_e;
            
            //double[] V2;
            //double[,] T_e;
            //double l1;
            //double m1;
            //double n1;
            //double l2;
            //double m2;
            //double n2;
            //double l3;
            //double m3;
            //double n3;

            
            //V2 = new double[3];
            //T_e = new double[6, 6];
            nGaussPoints = gp_d1 * gp_d2 * gp_d3; 
            ConsCartes = new double[nGaussPoints][,]; // TODOGer delete
            E = new double[nGaussPoints];
            ni = new double[nGaussPoints];
            //Cons = new double[6, 6];
            //Cons_T_e = new double[6, 6];
            for (int j = 0; j < nGaussPoints; j++)
            {
                double[] V3 = new double[3];
                double V3_norm;
                double[] V1 = new double[3];
                double V1_norm;
                
                
                E[j] = materialsAtGaussPoints[j].YoungModulus; // TODOGer delete
                ni[j] = materialsAtGaussPoints[j].PoissonRatio;
                //ConsCartes[j] = new double[6, 6];
                //for (int k = 0; k < 2; k++)
                //{ Cons[k, k] = E[j] / (1 - Math.Pow(ni[j], 2)); }
                //Cons[0, 1] = ni[j] * E[j] / (1 - Math.Pow(ni[j], 2));
                //Cons[1, 0] = ni[j] * E[j] / (1 - Math.Pow(ni[j], 2));
                //Cons[3, 3] = (1 - ni[j]) * (0.5) * E[j] / (1 - Math.Pow(ni[j], 2));
                //Cons[4, 4] = (1 - ni[j]) * (0.5) * E[j] / (1 - Math.Pow(ni[j], 2)); //Cons[4, 4] = (1 - ni[j]) * (0.41666666667) * E[j] / (1 - Math.Pow(ni[j], 2));
                //Cons[5, 5] = (1 - ni[j]) * (0.5) * E[j] / (1 - Math.Pow(ni[j], 2)); //Cons[5, 5] = (1 - ni[j]) * (0.41666666667) * E[j] / (1 - Math.Pow(ni[j], 2));
                ////for (int k = 0; k < 2; k++)
                ////{ Cons[4 + k, 4 + k] = (5 / 6) * (1 - ni[j]) * (0.5) * E[j] / (1 - Math.Pow(ni[j], 2)); }

                for (int k = 0; k < 3; k++)
                { V3[k] = 0; V1[k] = 0; /* V2[k] = 0; */ }

                for (int k = 0; k < 8; k++)
                {
                    for (int l = 0; l < 3; l++)
                    {
                        V3[l] += shapeFunctions[k][j] * oVn_i[k][l];
                        V1[l] += shapeFunctions[k][j] * oV1_i[k][l];
                    }
                }
                V3_norm = Math.Sqrt(V3[0] * V3[0] + V3[1] * V3[1] + V3[2] * V3[2]);
                V1_norm = Math.Sqrt(V1[0] * V1[0] + V1[1] * V1[1] + V1[2] * V1[2]);
                for (int l = 0; l < 3; l++)
                {
                    V3[l] = V3[l] / V3_norm;
                    V1[l] = V1[l] / V1_norm;
                }

                materialsAtGaussPoints[j].NormalVectorV3 = V3;
                materialsAtGaussPoints[j].TangentVectorV1 = V1;

                //V2[0] = V3[1] * V1[2] - V3[2] * V1[1];
                //V2[1] = V3[2] * V1[0] - V3[0] * V1[2];
                //V2[2] = V3[0] * V1[1] - V3[1] * V1[0];

                //l1 = V1[0];
                //m1 = V1[1];
                //n1 = V1[2];

                //l2 = V2[0];
                //m2 = V2[1];
                //n2 = V2[2];

                //l3 = V3[0];
                //m3 = V3[1];
                //n3 = V3[2];

                //for (int i = 0; i < 3; i++)
                //{
                //    T_e[0, i] = (V1[i] * V1[i]);
                //    T_e[1, i] = (V2[i] * V2[i]);
                //    T_e[2, i] = (V3[i] * V3[i]);

                //    T_e[3, i] = (2 * V1[i] * V2[i]);
                //    T_e[4, i] = (2 * V2[i] * V3[i]);
                //    T_e[5, i] = (2 * V3[i] * V1[i]);

                //    T_e[0, 3 + i] = (V1[i] * V1[1 + i - 3 * i * (i - 1) / 2]);
                //    T_e[1, 3 + i] = (V2[i] * V2[1 + i - 3 * i * (i - 1) / 2]);
                //    T_e[2, 3 + i] = (V3[i] * V3[1 + i - 3 * i * (i - 1) / 2]);

                //    T_e[3, 3 + i] = (V1[i] * V2[1 + i - 3 * i * (i - 1) / 2] + V2[i] * V1[1 + i - 3 * i * (i - 1) / 2]);
                //    T_e[4, 3 + i] = (V2[i] * V3[1 + i - 3 * i * (i - 1) / 2] + V3[i] * V2[1 + i - 3 * i * (i - 1) / 2]);
                //    T_e[5, 3 + i] = (V3[i] * V1[1 + i - 3 * i * (i - 1) / 2] + V1[i] * V3[1 + i - 3 * i * (i - 1) / 2]);
                //}

                // multiplication [Te']*[cons]*[Te];

                //for (int i = 0; i < 6; i++)
                //{
                //    for (int k = 0; k < 6; k++)
                //    {
                //        Cons_T_e[i, k] = 0;
                //        for (int l = 0; l < 6; l++)
                //        { Cons_T_e[i, k] += Cons[i, l] * T_e[l, k]; }
                //    }
                //}

                //for (int i = 0; i < 6; i++)
                //{
                //    for (int k = 0; k < 6; k++)
                //    {
                //        ConsCartes[j][i, k] = 0;
                //        for (int l = 0; l < 6; l++)
                //        { ConsCartes[j][i, k] += T_e[l, i] * Cons_T_e[l, k]; }
                //    }
                //}
            }

            //BL13 = GetBL13(shapeFunctionDerivatives, tUvec, J_0a);

            integrationCoefficient = new double[nGaussPoints];
            for (int j = 0; j < nGaussPoints; j++)
            { integrationCoefficient[j] += a_123g[j] * detJ_0[j]; }

        }

        

       
        //private double[,] SPK_circumflex;  //private double[][,] SPK_circumflex;

        //private void CalculateSPK()
        //{
        //    SPKvec = new double[nGaussPoints][];
        //    //SPK_circumflex = new double[9,9]; //SPK_circumflex = new double[nGaussPoints][,];
        //    for (int j = 0; j < nGaussPoints; j++)
        //    {
        //        SPKvec[j] = new double[6];
        //        //SPK_circumflex[j] = new double[9, 9];
        //        for (int l = 0; l < 6; l++)
        //        {
        //            SPKvec[j][l] = 0;
        //            for (int m = 0; m < 6; m++)
        //            {
        //                SPKvec[j][l] += ConsCartes[j][l, m] * GLvec[j][m];
        //            }

        //        }
        //        //for (int k = 0; k < 9; k++)
        //        //{
        //        //    for (int l = 0; l < 9; l++)
        //        //    { SPK_circumflex[j][k, l] = 0; }
        //        //}
        //        //for (int k = 0; k < 3; k++)
        //        //{
        //        //    for (int l = 0; l < 3; l++)
        //        //    {
        //        //        SPK_circumflex[j][3 * k + l, 3 * k + l] = SPKvec[j][l];
        //        //    }
        //        //    SPK_circumflex[j][3 * k, 3 * k + 1] = SPKvec[j][3];
        //        //    SPK_circumflex[j][3 * k, 3 * k + 2] = SPKvec[j][5];
        //        //    SPK_circumflex[j][3 * k + 1, 3 * k + 2] = SPKvec[j][4];

        //        //    SPK_circumflex[j][3 * k + 1, 3 * k] = SPKvec[j][3];
        //        //    SPK_circumflex[j][3 * k + 2, 3 * k] = SPKvec[j][5];
        //        //    SPK_circumflex[j][3 * k + 2, 3 * k + 1] = SPKvec[j][4];
        //        //}
        //    }
        //    //for (int k = 0; k < 9; k++)
        //    //{
        //    //    for (int l = 0; l < 9; l++)
        //    //    { SPK_circumflex[k, l] = 0; }  // PROSTHIKI RAM ektos loop nGausspoint
        //    //}
        //} 

        //private void UpdateSPK() //TODOGer Delete
        //{
        //    for (int j = 0; j < nGaussPoints; j++)
        //    {
        //        for (int l = 0; l < 6; l++)
        //        {
        //            SPKvec[j][l] = 0;
        //            for (int m = 0; m < 6; m++)
        //            {
        //                SPKvec[j][l] += ConsCartes[j][l, m] * GLvec[j][m];
        //            }

        //        }
        //        //for (int k = 0; k < 3; k++)
        //        //{
        //        //    for (int l = 0; l < 3; l++)
        //        //    {
        //        //        SPK_circumflex[j][3 * k + l, 3 * k + l] = SPKvec[j][l];
        //        //    }
        //        //    SPK_circumflex[j][3 * k, 3 * k + 1] = SPKvec[j][3];
        //        //    SPK_circumflex[j][3 * k, 3 * k + 2] = SPKvec[j][5];
        //        //    SPK_circumflex[j][3 * k + 1, 3 * k + 2] = SPKvec[j][4];

        //        //    SPK_circumflex[j][3 * k + 1, 3 * k] = SPKvec[j][3];
        //        //    SPK_circumflex[j][3 * k + 2, 3 * k] = SPKvec[j][5];
        //        //    SPK_circumflex[j][3 * k + 2, 3 * k + 1] = SPKvec[j][4];
        //        //}

        //    }
        //}

        //private double[][,] ck;// 1 ana komvo kai ana gauss Point dld [GP][8komvoi,diastash9]
        private double [][,] CalculateCk(double[] [,] J_0a,double [][] tU )
        {
            //initialize
            nGaussPoints = gp_d1 * gp_d2 * gp_d3;
            double[][,] ck;
            ck = new double[nGaussPoints][,];
            for (int j = 0; j < nGaussPoints; j++)
            {
                ck[j] = new double[8, 9];
            }
            //tupoi
            for (int j = 0; j < nGaussPoints; j++)
            {
                for (int m = 0; m < 8; m++)
                {
                    for (int k = 0; k < 3; k++)
                    {
                        for (int l = 0; l < 3; l++)
                        {
                            ck[j][m, 3 * k + l] = J_0a[j][l, 2 * m + 1] * tU[m][3 + k];
                        }
                    }
                }
            }
            return ck;
        }

        //private void UpdateCk()
        //{
        //    for (int j = 0; j < nGaussPoints; j++)
        //    {
        //        for (int m = 0; m < 8; m++)
        //        {
        //            for (int k = 0; k < 3; k++)
        //            {
        //                for (int l = 0; l < 3; l++)
        //                {
        //                    ck[j][m, 3 * k + l] = J_0a[j][l, 2 * m + 1] * tU[m][3 + k];
        //                }
        //            }
        //        }
        //    }
        //}



        //private void UpdatePartiallyPrecalculatedVariables__forStrains(Element element, double[][] tx_i)
        //{
        //    //this.UpdateJ_1b();
        //    //this.UpdateJ_1();
        //    //this.UpdateDefGradTr();
        //    //this.UpdateGL();
        //    //this.UpdateGLvec();
        //    this.CalculateStrains(element, tx_i);
        //}

        private void CalculateStrains(IElement element, double[][] tx_i)
        {
            // prosthiki logw J_0inv oxi global
            // kai anagkastika kai J_0a
            (double[][] gausscoordinates, double[][] shapeFunctions, double[][] shapeFunctionDerivatives, double[] a_123g) =
                interpolation.GetShapeFunctions(gp_d1, gp_d2, gp_d3);
            (double[][,] ll1, double[][,] J_0a) = JacobianShell8Calculations.Getll1AndJ_0a(
                nGaussPoints, tk, gausscoordinates, shapeFunctions, shapeFunctionDerivatives);


            (double[][,] J_0inv, double[] detJ_0) =
                JacobianShell8Calculations.GetJ_0invAndDetJ_0(J_0a, element.INodes, oVn_i, nGaussPoints);

            double[,] J_1b;
            double[][,] J_1;
            double[][,] DefGradTr;
            double[][,] GL;

            J_1b = new double[16, 3];
            J_1 = new double[nGaussPoints][,];
            DefGradTr = new double[nGaussPoints][,];
            GL = new double[nGaussPoints][,];
            for (int j = 0; j < nGaussPoints; j++)
            {
                J_1[j] = new double[3, 3];
                DefGradTr[j] = new double[3, 3];
                GL[j] = new double[3, 3];
            }

            GLvec = new double[nGaussPoints][]; // xehwristh perioxh 

            for (int j = 0; j < 8; j++)
            {
                J_1b[2 * j, 0] = tx_i[j][0];
                J_1b[2 * j + 1, 0] = tU[j][3];
                J_1b[2 * j, 1] = tx_i[j][1];
                J_1b[2 * j + 1, 1] = tU[j][4];
                J_1b[2 * j, 2] = tx_i[j][2];
                J_1b[2 * j + 1, 2] = tU[j][5];
            }

            for (int j = 0; j < nGaussPoints; j++)
            {
                for (int k = 0; k < 3; k++)
                {
                    for (int l = 0; l < 3; l++)
                    {
                        J_1[j][k, l] = 0;
                        for (int m = 0; m < 16; m++)
                        {
                            J_1[j][k, l] += J_0a[j][k, m] * J_1b[m, l];
                        }

                    }

                }

                for (int k = 0; k < 3; k++)
                {
                    for (int l = 0; l < 3; l++)
                    {
                        DefGradTr[j][k, l] = 0;
                        for (int m = 0; m < 3; m++)
                        {
                            DefGradTr[j][k, l] += J_0inv[j][k, m] * J_1[j][m, l];
                        }

                    }

                }

                for (int k = 0; k < 3; k++)
                {
                    for (int l = 0; l < 3; l++)
                    {
                        GL[j][k, l] = 0;
                        for (int m = 0; m < 3; m++)
                        {
                            GL[j][k, l] += DefGradTr[j][k, m] * DefGradTr[j][l, m];
                        }

                    }

                }
                for (int k = 0; k < 3; k++)
                {
                    GL[j][k, k] = GL[j][k, k] - 1;
                }
                for (int k = 0; k < 3; k++)
                {
                    for (int l = 0; l < 3; l++)
                    { GL[j][k, l] = 0.5 * GL[j][k, l]; }
                }


                GLvec[j] = new double[6];
                for (int k = 0; k < 3; k++)
                { GLvec[j][k] = GL[j][k, k]; }
                GLvec[j][3] = 2 * GL[j][0, 1];
                GLvec[j][4] = 2 * GL[j][1, 2];
                GLvec[j][5] = 2 * GL[j][2, 0];

            }
        }

        //private void UpdatePartiallyPrecalculatedVariables_andforForces()
        //{
        //    this.Updatell2();
        //    this.Updatel_circumflex();
        //    this.UpdateBL11b();
        //    this.UpdateBL11();
        //    this.UpdateBL13();

        //}

        //private void UpdatePartiallyPrecalculatedVariables_andforStiffnessMatrix()
        //{
        //    this.UpdateCk();
        //    //mporei na periexei kai SPK_circumflex apo to aplo ean xreiasthei
        //}



        //private double[][] kck;// 1 ana komvo kai (ana Gauss point+1 gia to athroiskma) [GP][8 vathmoi komvoi]
                               // to initialize tou einai comment out parapanw

        //private double[,] BNL;  //private double[][,] BNL; //PROSTHIKI RAM
        //private double[][,] KNL;

        //private double[][,] KL;
        //private double[][,] BL1_2;
        //private double[][,] BL1;
        //private double[][,] BL0;
        //private double[][,] BL;

        //private double[][,] ConsBL;
        //private double[][,] S_BNL;

        //private double[][,] BL01plus1_2;
        //private double[][] BL01plus1_2tSPKvec;

        //private double[,] Kt = new double[40, 40];

        //private double[][] Fxk;

        private void InitializeFandKmatrixes()
        {
            //BNL = new double[9, 40];  //BNL = new double[nGaussPoints][,];
            //BL1_2 = new double[nGaussPoints][,];
            //BL1 = new double[nGaussPoints][,];
            //BL0 = new double[nGaussPoints][,];
            //BL = new double[nGaussPoints][,];

            //kck = new double[nGaussPoints + 1][];
            //KL = new double[nGaussPoints + 1][,];
            //KNL = new double[nGaussPoints + 1][,];

            //ConsBL = new double[nGaussPoints][,];
            //S_BNL = new double[nGaussPoints][,];

            //kck = new double[nGaussPoints + 1][];

            //BL01plus1_2 = new double[nGaussPoints][,];
            //BL01plus1_2tSPKvec = new double[nGaussPoints][];

            //Kt = new double[40, 40];

            //Fxk = new double[nGaussPoints + 1][];

            for (int j = 0; j < nGaussPoints; j++)
            {
                //BNL[j] = new double[9, 40]; // PROSTHIKI RAM 
                //BL1_2[j] = new double[6, 9];
                //BL1[j] = new double[6, 40];
                //BL0[j] =new double[6, 40];
                //BL[j] = new double[6, 40];

                //ConsBL[j] = new double[6, 40];
                //S_BNL[j] = new double[9, 40];

                //kck[j] = new double[8];

                //BL01plus1_2[j] = new double[6, 9];
                //BL01plus1_2tSPKvec[j] = new double[9];
            }

            for (int j = 0; j < nGaussPoints + 1; j++)
            {
                //kck[j] = new double[8];
                //KL[j] = new double[40, 40];
                //KNL[j] = new double[40, 40];

                //Fxk[j] = new double[40];
            }
        }

        private double [] UpdateForces(IElement element)
        {

            //prosthikes gia ll1 entos methodwn mono
            (double[][] gausscoordinates, double[][] shapeFunctions, double[][] shapeFunctionDerivatives, double[] a_123g) =
                interpolation.GetShapeFunctions(gp_d1, gp_d2, gp_d3);
            (double[][,] ll1, double[][,] J_0a) = JacobianShell8Calculations.Getll1AndJ_0a(
                nGaussPoints, tk, gausscoordinates, shapeFunctions, shapeFunctionDerivatives);
            //prosthikes gia ll1 entos methodwn mono

            //prosthikes gia BL11a entos methodwn
            (double[][,] J_0inv, double[] detJ_0) =
                JacobianShell8Calculations.GetJ_0invAndDetJ_0(J_0a, element.INodes, oVn_i, nGaussPoints);
            double[][,] BL11a;
            BL11a = GetBL11a(J_0inv);
            double[][,] BL12;
            BL12 = GetBL12(J_0inv);
            double[][,] BL13;
            BL13 = GetBL13(shapeFunctionDerivatives, tUvec, J_0a);

            //prosthiki declare BL klp entos upd force
            double[][,] BL;
            BL = new double[nGaussPoints][,];
            double[][,] BL01plus1_2;
            BL01plus1_2 = new double[nGaussPoints][,];
            for (int j = 0; j < nGaussPoints; j++)
            {
                BL[j] = new double[6, 40];
                BL01plus1_2[j] = new double[6, 9];
            }


            double[][] Fxk= new double [nGaussPoints+1] [];
            for (int j = 0; j<nGaussPoints + 1; j++)
            {
                Fxk[j] = new double[40];
            }

            // PROSTHIKI apo declare ektos klashs
             double[,] ll2;
            ll2 = new double[24, 3];                 
            for (int j = 0; j < 8; j++)
            {
                for (int k = 0; k < 3; k++)
                {
                    ll2[3 * j + 0, k] = tU[j][k];
                    ll2[3 * j + 1, k] = tU[j][3 + k];
                    ll2[3 * j + 2, k] = oVn_i[j][k];
                }
            }

            double[][,] l_circumflex;
            l_circumflex = new double[nGaussPoints][,];
            for (int j = 0; j < nGaussPoints; j++)
            { l_circumflex[j] = new double[3, 3]; }
            for (int j = 0; j < nGaussPoints; j++)
            {
                for (int k = 0; k < 3; k++)
                {
                    for (int l = 0; l < 3; l++)
                    {
                        l_circumflex[j][k, l] = 0;
                        for (int m = 0; m < 24; m++)
                        {
                            l_circumflex[j][k, l] += ll1[j][k, m] * ll2[m, l];
                        }

                    }

                }

            }

            double[][,] BL11b;
            BL11b = new double[nGaussPoints][,];
            for (int j = 0; j < nGaussPoints; j++)
            { BL11b[j] = new double[9, 9]; }
            for (int j = 0; j < nGaussPoints; j++)
            {
                for (int k = 0; k < 9; k++)
                {
                    for (int l = 0; l < 9; l++)
                    { BL11b[j][k, l] = 0; }
                }


                for (int k = 0; k < 3; k++)
                {
                    for (int l = 0; l < 3; l++)
                    {
                        for (int m = 0; m < 3; m++)
                        { BL11b[j][3 * k + l, 3 * k + m] = l_circumflex[j][l, m]; }
                    }
                }
            }

            double[][,] BL11;
            BL11 = new double[nGaussPoints][,];
            for (int j = 0; j < nGaussPoints; j++)
            { BL11[j] = new double[6, 9]; }
            for (int j = 0; j < nGaussPoints; j++)
            {
                for (int k = 0; k < 6; k++)
                {
                    for (int l = 0; l < 9; l++)
                    {
                        BL11[j][k, l] = 0;
                        for (int m = 0; m < 9; m++)
                        {
                            BL11[j][k, l] += BL11a[j][k, m] * BL11b[j][m, l];
                        }
                    }
                }
            }

            
            double[][,] BL1_2;
            BL1_2 = new double[nGaussPoints][,];
            for (int j = 0; j < nGaussPoints; j++)
            {
                BL1_2[j] = new double[6, 9];
            }
            for (int j = 0; j < nGaussPoints; j++)
            {
                for (int k = 0; k < 6; k++)
                {
                    for (int l = 0; l < 9; l++)
                    {
                        BL1_2[j][k, l] = 0;
                        for (int m = 0; m < 9; m++)
                        {
                            BL1_2[j][k, l] += BL11[j][k, m] * BL12[j][m, l]; //TODO BL11 keno kai BL12 null thelei getbl12 edw kai parakatw
                        }                                                   //vriskomaste sto calculate Kmatrices eprepe na trexei to calculate BL11 prwta

                    }

                }
                //for (int k = 0; k < 6; k++)
                //{
                //    for (int l = 0; l < 40; l++)
                //    {
                //        BL1[j][k, l] = 0;
                //        BL0[j][k, l] = 0;
                //        for (int m = 0; m < 9; m++) //panw apo to for BLx=BL1_2+BL11 kai mesa sto for BL=BLx*BL13
                //        {
                //            BL1[j][k, l] += BL1_2[j][k, m] * BL13[j][m, l];
                //            BL0[j][k, l] += BL11[j][k, m]* BL13[j][m, l];
                //        }
                //        BL[j][k, l] = BL0[j][k, l] + BL1[j][k, l];
                //    }
                //}

                for (int k = 0; k < 6; k++)
                {
                    for (int l = 0; l < 9; l++)
                    {
                        BL01plus1_2[j][k, l] = BL1_2[j][k, l] + BL11a[j][k, l];
                    }
                }

                for (int k = 0; k < 6; k++)
                {
                    for (int l = 0; l < 40; l++)
                    {
                        BL[j][k, l] = 0;
                        for (int m = 0; m < 9; m++)
                        {
                            BL[j][k, l] += BL01plus1_2[j][k, m] * BL13[j][m, l];
                        }
                    }
                }
            }

            //mprfwsi drasewn   
            for (int j = 0; j < nGaussPoints; j++)
            {
                for (int k = 0; k < 40; k++)
                {
                    Fxk[j][k] = 0;
                    for (int m = 0; m < 6; m++)
                    { Fxk[j][k] += BL[j][m, k] * materialsAtGaussPoints[j].Stresses[m]; }
                }
            }
            for (int k = 0; k < 40; k++)
            {
                Fxk[nGaussPoints][k] = 0;
                for (int j = 0; j < nGaussPoints; j++)
                { Fxk[nGaussPoints][k] += integrationCoefficient[j] * Fxk[j][k]; }
            }

            return Fxk[nGaussPoints];
        }

        private double [,] UpdateKmatrices(IElement element)
        {
            double[,] Kt = new double[40, 40];
            // PROSTHIKI RAM apo osa declared ektos methodou
            double[,] SPK_circumflex;  //private double[][,] SPK_circumflex;
            double[,] BNL;
            double[,] S_BNL;
            double[,] ConsBL;
            double[][] kck;// 1 ana komvo kai (ana Gauss point+1 gia to athroiskma) [GP][8 vathmoi komvoi]
            double[][,] KNL;
            double[][,] KL;
            double[][] BL01plus1_2tSPKvec;



            (double[][] gausscoordinates, double[][] shapeFunctions, double[][] shapeFunctionDerivatives, double[] a_123g) =
                interpolation.GetShapeFunctions(gp_d1, gp_d2, gp_d3);
            (double[][,] ll1, double[][,] J_0a) = JacobianShell8Calculations.Getll1AndJ_0a(
                nGaussPoints, tk, gausscoordinates, shapeFunctions, shapeFunctionDerivatives);
            (double[][,] J_0inv, double[] detJ_0) =
                JacobianShell8Calculations.GetJ_0invAndDetJ_0(J_0a, element.INodes, oVn_i, nGaussPoints);
            double[][,] BNL1;
            BNL1 = GetBNL1(J_0inv);
            double[][,] BL13;
            BL13 = GetBL13(shapeFunctionDerivatives, tUvec, J_0a); //logw pou de tha kratietai apo Ton upologismo tou stis forces



            double[][,] ck;
            ck = CalculateCk(J_0a, tU);


            //PERIOXH APO FORCES LOGW BL kai BL01plus1_2 den tha kratiountai 
            // kai anagkastika kai ta B:L11a kai BL12 prwta
            double[][,] BL11a;
            BL11a = GetBL11a(J_0inv);
            double[][,] BL12;
            BL12 = GetBL12(J_0inv);
            // kai ana gkastika ta declare kai ta initializations
            double[][,] BL;
            BL = new double[nGaussPoints][,];
            double[][,] BL01plus1_2;
            BL01plus1_2 = new double[nGaussPoints][,];
            for (int j = 0; j < nGaussPoints; j++)
            {
                BL[j] = new double[6, 40];
                BL01plus1_2[j] = new double[6, 9];
            }
            //
            double[,] ll2;
            ll2 = new double[24, 3];
            for (int j = 0; j < 8; j++)
            {
                for (int k = 0; k < 3; k++)
                {
                    ll2[3 * j + 0, k] = tU[j][k];
                    ll2[3 * j + 1, k] = tU[j][3 + k];
                    ll2[3 * j + 2, k] = oVn_i[j][k];
                }
            }

            double[][,] l_circumflex;
            l_circumflex = new double[nGaussPoints][,];
            for (int j = 0; j < nGaussPoints; j++)
            { l_circumflex[j] = new double[3, 3]; }
            for (int j = 0; j < nGaussPoints; j++)
            {
                for (int k = 0; k < 3; k++)
                {
                    for (int l = 0; l < 3; l++)
                    {
                        l_circumflex[j][k, l] = 0;
                        for (int m = 0; m < 24; m++)
                        {
                            l_circumflex[j][k, l] += ll1[j][k, m] * ll2[m, l];
                        }

                    }

                }

            }

            double[][,] BL11b;
            BL11b = new double[nGaussPoints][,];
            for (int j = 0; j < nGaussPoints; j++)
            { BL11b[j] = new double[9, 9]; }
            for (int j = 0; j < nGaussPoints; j++)
            {
                for (int k = 0; k < 9; k++)
                {
                    for (int l = 0; l < 9; l++)
                    { BL11b[j][k, l] = 0; }
                }


                for (int k = 0; k < 3; k++)
                {
                    for (int l = 0; l < 3; l++)
                    {
                        for (int m = 0; m < 3; m++)
                        { BL11b[j][3 * k + l, 3 * k + m] = l_circumflex[j][l, m]; }
                    }
                }
            }

            double[][,] BL11;
            BL11 = new double[nGaussPoints][,];
            for (int j = 0; j < nGaussPoints; j++)
            { BL11[j] = new double[6, 9]; }
            for (int j = 0; j < nGaussPoints; j++)
            {
                for (int k = 0; k < 6; k++)
                {
                    for (int l = 0; l < 9; l++)
                    {
                        BL11[j][k, l] = 0;
                        for (int m = 0; m < 9; m++)
                        {
                            BL11[j][k, l] += BL11a[j][k, m] * BL11b[j][m, l];
                        }
                    }
                }
            }
            double[][,] BL1_2;
            BL1_2 = new double[nGaussPoints][,];
            for (int j = 0; j < nGaussPoints; j++)
            {
                BL1_2[j] = new double[6, 9];
            }
            for (int j = 0; j < nGaussPoints; j++)
            {
                for (int k = 0; k < 6; k++)
                {
                    for (int l = 0; l < 9; l++)
                    {
                        BL1_2[j][k, l] = 0;
                        for (int m = 0; m < 9; m++)
                        {
                            BL1_2[j][k, l] += BL11[j][k, m] * BL12[j][m, l]; //TODO BL11 keno kai BL12 null thelei getbl12 edw kai parakatw
                        }                                                   //vriskomaste sto calculate Kmatrices eprepe na trexei to calculate BL11 prwta

                    }

                }
                for (int k = 0; k < 6; k++)
                {
                    for (int l = 0; l < 9; l++)
                    {
                        BL01plus1_2[j][k, l] = BL1_2[j][k, l] + BL11a[j][k, l];
                    }
                }

                for (int k = 0; k < 6; k++)
                {
                    for (int l = 0; l < 40; l++)
                    {
                        BL[j][k, l] = 0;
                        for (int m = 0; m < 9; m++)
                        {
                            BL[j][k, l] += BL01plus1_2[j][k, m] * BL13[j][m, l];
                        }
                    }
                }
            }
            // PERIOXH APO FORCES LOGW BL kai BL01plu1_2 den tha kratiountai 


            // PROSTHIKI RAM 2 apo osa ginontai initialized entos methodou
            SPK_circumflex = new double[9, 9];
            BNL = new double[9, 40];
            S_BNL = new double[9, 40];
            ConsBL = new double[6, 40];
            KL = new double[nGaussPoints + 1][,];
            KNL = new double[nGaussPoints + 1][,];
            kck = new double[nGaussPoints + 1][];
            BL01plus1_2tSPKvec = new double[nGaussPoints][];
            for (int j = 0; j < nGaussPoints + 1; j++)
            {
                kck[j] = new double[8];
                KL[j] = new double[40, 40];
                KNL[j] = new double[40, 40];
            }
            for (int j = 0; j < nGaussPoints; j++)
            {
                BL01plus1_2tSPKvec[j] = new double[9];
            }


            for (int k = 0; k < 9; k++)
            {
                for (int l = 0; l < 9; l++)
                { SPK_circumflex[k, l] = 0; }  // PROSTHIKI RAM ektos loop nGausspoint
            }
            for (int j = 0; j < nGaussPoints; j++)
            {
                for (int k = 0; k < 3; k++) // PROSTHIKI RAM den kratietai to 
                {
                    for (int l = 0; l < 3; l++)
                    {
                        SPK_circumflex[3 * k + l, 3 * k + l] = materialsAtGaussPoints[j].Stresses[l];
                    }
                    SPK_circumflex[3 * k, 3 * k + 1] = materialsAtGaussPoints[j].Stresses[3];
                    SPK_circumflex[3 * k, 3 * k + 2] = materialsAtGaussPoints[j].Stresses[5];
                    SPK_circumflex[3 * k + 1, 3 * k + 2] = materialsAtGaussPoints[j].Stresses[4];

                    SPK_circumflex[3 * k + 1, 3 * k] = materialsAtGaussPoints[j].Stresses[3];
                    SPK_circumflex[3 * k + 2, 3 * k] = materialsAtGaussPoints[j].Stresses[5];
                    SPK_circumflex[3 * k + 2, 3 * k + 1] = materialsAtGaussPoints[j].Stresses[4];
                }

                for (int k = 0; k < 9; k++)
                {
                    for (int l = 0; l < 40; l++)
                    {
                        BNL[k, l] = 0;
                        for (int m = 0; m < 9; m++)
                        {
                            BNL[k, l] += BNL1[j][k, m] * BL13[j][m, l];
                        }

                    }
                }

                for (int k = 0; k < 9; k++)
                {
                    BL01plus1_2tSPKvec[j][k] = 0;
                    for (int m = 0; m < 6; m++)
                    {
                        BL01plus1_2tSPKvec[j][k] += BL01plus1_2[j][m, k] * materialsAtGaussPoints[j].Stresses[m];
                    }
                }

                for (int k = 0; k < 8; k++)
                {
                    kck[j][k] = 0;
                    for (int m = 0; m < 9; m++)
                    {
                        kck[j][k] += ck[j][k, m] * BL01plus1_2tSPKvec[j][m];
                    }
                }

                // porsthetoume kai to kck ws extra(den prokuptei apo ta comment out

                for (int k = 0; k < 6; k++)
                {
                    for (int l = 0; l < 40; l++)
                    {
                        ConsBL[k, l] = 0;
                        for (int m = 0; m < 6; m++)
                        {
                            ConsBL[k, l] += ConsCartes[j][k, m] * BL[j][m, l];
                        }
                    }
                }

                for (int k = 0; k < 9; k++)
                {
                    for (int l = 0; l < 40; l++)
                    {
                        S_BNL[k, l] = 0;
                        for (int m = 0; m < 9; m++)
                        {
                            S_BNL[k, l] += SPK_circumflex[k, m] * BNL[m, l];
                        }
                    }
                }

                for (int k = 0; k < 40; k++)
                {
                    for (int l = 0; l < 40; l++)
                    {
                        KNL[j][k, l] = 0;
                        for (int m = 0; m < 9; m++)
                        {
                            KNL[j][k, l] += BNL[m, k] * S_BNL[m, l];
                        }
                    }
                }

                for (int k = 0; k < 40; k++)
                {
                    for (int l = 0; l < 40; l++)
                    {
                        KL[j][k, l] = 0;
                        for (int m = 0; m < 6; m++)
                        {
                            KL[j][k, l] += BL[j][m, k] * ConsBL[m, l];
                        }
                    }
                }
            }

            // morfwsi telikou mhtrwou
            for (int k = 0; k < 40; k++)
            {
                for (int l = 0; l < 40; l++)
                { Kt[k, l] = 0; }
            }
            for (int j = 0; j < nGaussPoints; j++)
            {

                for (int k = 0; k < 40; k++)
                {
                    for (int l = 0; l < 40; l++)
                    { Kt[k, l] += integrationCoefficient[j] * (KL[j][k, l] + KNL[j][k, l]); }
                }

                for (int l = 0; l < 8; l++)
                {
                    Kt[5 * l + 3, 5 * l + 3] += integrationCoefficient[j] * kck[j][l];
                    Kt[5 * l + 4, 5 * l + 4] += integrationCoefficient[j] * kck[j][l];
                }
            }

            return Kt;
        }

        // ANANEWSH thw thesis tou stoixeiou-----------------------------------------
        // voithitikes metavlhtes gia upologismo strofhs-----------------------------
        private double[] ak_total = new double[8];
        private double[] bk_total = new double[8];

        // metavlhtes gia anafora stis strofes kai voithitikoi pinakes
        //private double ak;
        //private double bk;
        //private double gk1;
        //private double[,] Q = new double[3, 3];
        //private double[,] Q2 = new double[3, 3];
        //private double[] tdtVn = new double[3];
        //private double tV1norm;

        private void UpdateCoordinateData(double[] localdisplacements, IList<INode> elementNodes, out double[][] tx_i)
        {
            double[][] ox_i = new double[8][]; // this should be allocated in the constructor
            for (int j = 0; j < 8; j++)
            {
                ox_i[j] = new double[] { elementNodes[j].X, elementNodes[j].Y, elementNodes[j].Z, };
            }

            //PROSHIKI ram
            double ak;
            double bk;
            tx_i = new double[8][];
            for (int k = 0; k < 8; k++)
            {
                tx_i[k] = new double[3];
                for (int l = 0; l < 3; l++)
                {
                    tx_i[k][l] = ox_i[k][l] + localdisplacements[5 * k + l];
                }
                //update twn tU kai tUvec 
                tU[k][0] = localdisplacements[5 * k + 0];
                tU[k][1] = localdisplacements[5 * k + 1];
                tU[k][2] = localdisplacements[5 * k + 2];
                ak = localdisplacements[5 * k + 3] - ak_total[k];
                ak_total[k] = localdisplacements[5 * k + 3];
                bk = localdisplacements[5 * k + 4] - bk_total[k];
                bk_total[k] = localdisplacements[5 * k + 4];
                this.RotateNodalDirectionVectors(ak, bk, k);
                // update twn tU kai tUvec ews edw         
            }
            // shmeio print dedomenwn gia debug
            ////PrintUtilities.ConvertAndWriteToFileVector(tU, @"C:\Users\turbo-x\Desktop\cohesive_check_MSOLVE_2\paradeigma_apo_arxika_swsta_shell_orthi_gia_check_tou_neou_class\orthi\CopyApoTaShellNewLoadCaseArgurhs\tU_SHELL_local_msolve1.txt");
            ////PrintUtilities.ConvertAndWriteToFileVector(tUvec, @"C:\Users\turbo-x\Desktop\cohesive_check_MSOLVE_2\paradeigma_apo_arxika_swsta_shell_orthi_gia_check_tou_neou_class\orthi\CopyApoTaShellNewLoadCaseArgurhs\tUvec_SHELL_local_msolve1.txt");
        }

        // voithitikes metavlhtes gia thn peristrofh
        //private double[] tdtVn = new double[3];
        //private double[] tdtV1 = new double[3];
        //private double[] tdtV2 = new double[3];
        //private double theta;
        //private double[] theta_vec = new double[3];
        //private double[,] s_k = new double[3, 3];
        private void RotateNodalDirectionVectors(double ak, double bk, int n_vector)
        {
            //PROSTHIKI gia ram
            double gk1;
            double[,] Q = new double[3, 3];
            double[,] Q2 = new double[3, 3];
            double[] tdtVn = new double[3];
            double[] tdtV1 = new double[3];
            double[] tdtV2 = new double[3];
            double theta;
            double[] theta_vec = new double[3];
            double[,] s_k = new double[3, 3];


            for (int j = 0; j < 3; j++)
            {
                theta_vec[j] = ak * tUvec[n_vector][j] + bk * tUvec[n_vector][3 + j];
            }
            theta = Math.Sqrt((theta_vec[0] * theta_vec[0]) + (theta_vec[1] * theta_vec[1]) + (theta_vec[2] * theta_vec[2]));
            if (theta > 0)
            {
                s_k[0, 1] = -theta_vec[2];
                s_k[0, 2] = theta_vec[1];
                s_k[1, 0] = theta_vec[2];
                s_k[1, 2] = -theta_vec[0];
                s_k[2, 0] = -theta_vec[1];
                s_k[2, 1] = theta_vec[0];

                for (int j = 0; j < 3; j++)
                {
                    for (int m = 0; m < 3; m++)
                    {
                        Q[j, m] = (Math.Sin(theta) / theta) * s_k[j, m];
                    }
                }

                for (int m = 0; m < 3; m++)
                {
                    Q[m, m] += 1;
                }
                gk1 = 0.5 * ((Math.Sin(0.5 * theta) / (0.5 * theta)) * (Math.Sin(0.5 * theta) / (0.5 * theta)));
                for (int j = 0; j < 3; j++)
                {
                    for (int m = 0; m < 3; m++)
                    {
                        Q2[j, m] = 0;
                        for (int n = 0; n < 3; n++)
                        { Q2[j, m] += gk1 * s_k[j, n] * s_k[n, m]; }
                    }
                }
                for (int j = 0; j < 3; j++)
                {
                    for (int m = 0; m < 3; m++)
                    {
                        Q[j, m] += Q2[j, m];
                    }
                }
                //
                for (int j = 0; j < 3; j++)
                {
                    tdtVn[j] = 0;
                    for (int m = 0; m < 3; m++)
                    {
                        tdtVn[j] += Q[j, m] * tU[n_vector][3 + m];
                    }
                }

                for (int j = 0; j < 3; j++)
                {
                    tU[n_vector][3 + j] = tdtVn[j];
                }
                //
                for (int j = 0; j < 3; j++)
                {
                    tdtV1[j] = 0;
                    for (int m = 0; m < 3; m++)
                    {
                        tdtV1[j] += Q[j, m] * tUvec[n_vector][m];
                    }
                }

                for (int j = 0; j < 3; j++)
                {
                    tUvec[n_vector][j] = tdtV1[j];
                }
                //
                for (int j = 0; j < 3; j++)
                {
                    tdtV2[j] = 0;
                    for (int m = 0; m < 3; m++)
                    {
                        tdtV2[j] += Q[j, m] * tUvec[n_vector][3 + m];
                    }
                }

                for (int j = 0; j < 3; j++)
                {
                    tUvec[n_vector][3 + j] = tdtV2[j];
                }
            }
        }

        // aparaithta tou IStructuralFiniteElement

        public int ID
        {
            get { return 12; }
        }

        public ElementDimensions ElementDimensions
        {
            get { return ElementDimensions.ThreeD; }
        }

        public virtual IList<IList<DOFType>> GetElementDOFTypes(IElement element)
        {
            return dofTypes;
        }

        public IElementDOFEnumerator DOFEnumerator
        {
            get { return dofEnumerator; }
            set { dofEnumerator = value; }
        }

        //aparaithta tou IstructuralElement gia to material
        public void ClearMaterialState()
        {
            foreach (IIsotropicContinuumMaterial3D m in materialsAtGaussPoints) m.ClearState();
        }

        public void SaveMaterialState()
        {
            foreach (IIsotropicContinuumMaterial3D m in materialsAtGaussPoints) m.SaveState();
        }

        public void ClearMaterialStresses()
        {
            foreach (IIsotropicContinuumMaterial3D m in materialsAtGaussPoints) m.ClearStresses();
        }

        public bool MaterialModified
        {
            get
            {
                foreach (IIsotropicContinuumMaterial3D material in materialsAtGaussPoints)
                    if (material.Modified) return true;
                return false;
            }
        }

        public void ResetMaterialModified()
        {
            foreach (IIsotropicContinuumMaterial3D material in materialsAtGaussPoints) material.ResetModified();
        }

        public Tuple<double[], double[]> CalculateStresses(Element element, double[] localTotalDisplacements, double[] localdDisplacements)
        {
            this.UpdateCoordinateData(localTotalDisplacements, element.INodes, out double[][] tx_i);
            this.CalculateStrains(element, tx_i);
            for (int npoint = 0; npoint < materialsAtGaussPoints.Length; npoint++)
            {
                materialsAtGaussPoints[npoint].UpdateMaterial(GLvec[npoint]);
            }
            //this.UpdatePartiallyPrecalculatedVariables__forStrains( element, tx_i);
            //this.UpdateSPK(); //mporei na lamvanetai apo uliko nme materialsAtGPs.Stresses // mporei na xwristhei afto se spkvec kai ta upoloiopa pou einai gia KMatrices (SPK_circumflex)

            return new Tuple<double[], double[]>(new double[123], materialsAtGaussPoints[materialsAtGaussPoints.Length - 1].Stresses);
        }

        //aparaithta tou Istructural gia th dunamiki analusi
        public double[] CalculateAccelerationForces(Element element, IList<MassAccelerationLoad> loads)
        {
            return new double[123];
        }

        public virtual IMatrix2D MassMatrix(IElement element)
        {
            throw new NotImplementedException();
        }

        public virtual IMatrix2D DampingMatrix(IElement element)
        {
            throw new NotImplementedException();
        }

        // forces tha exei ena if me initial geometric data
        // kai to else me strofi kai meta olous tous upologismous.

        // mporei kai me ena if ean einai to trito dianusma mhdeniko. na mhn ektelountai kapoioi upologismoi.

        //implementation of basic methods

        public double[] CalculateForces(Element element, double[] localTotalDisplacements, double[] localdDisplacements)
        {
            double[] Fxk;
           
            //this.UpdatePartiallyPrecalculatedVariables_andforForces();

            // TODOGer metafora sto spkvec to material.stresses h xrhsh kateftheian stresses ** ti ginetai me to circumflex
            Fxk=this.UpdateForces(element);

            //PrintUtilities.WriteToFileVector(Fxk, @"C:\Users\turbo-x\Desktop\notes_elegxoi\MSOLVE_output\Fxk_shell.txt");
            return Fxk;
        }

        private int endeixiStiffness = 1;
        public virtual IMatrix2D StiffnessMatrix(IElement element)
        {
            double[,] Kt = new double[40, 40];
            if (endeixiStiffness == 1)
            {
                //this.CalculateGaussCoordinatesShapefunctionDataAndll1J0_a();
                //GetInitialGeometricData(element);
                //this.CalculateCompletelyPrecalculatedVariables();
                //this.CalculateCons();
                this.CalculateInitialConfigurationData(element, out double[][] tx_i);


                //this.InitializeAndCalculateOriginalValuesForPartiallyPrecalculatedVariables(); //<-- periexei to calculate strains                
                this.CalculateStrains(element, tx_i);

                //this.CalculateCk();
                //this.CalculateSPK(); // ousiastika kanei th douleia tou ulikou
                //edw to hexa pairnei tis Spkvec[npoint][j] = materialsAtGaussPoints[npoint].Stresses[j];
                //TODOGer tha to kaname etsi alla pleon tha xrhs mono materials.stresses kateftheian.
                //Gia arxikes strains: (TODOGer den exei ginei sto hexa8)
                for (int npoint = 0; npoint < materialsAtGaussPoints.Length; npoint++)
                {
                    materialsAtGaussPoints[npoint].UpdateMaterial(GLvec[npoint]);
                }

                this.InitializeFandKmatrixes();

                //this.Calculatell2();
                //this.Calculatel_circumflex();
                //this.CalculateBL11b();
                //this.CalculateBL11();
                //this.CalculateBL13();
                this.UpdateForces(element); // dioti periexei kai mhtrwa BL pou einai aparaithta gia to stifness matrix


                Kt=this.UpdateKmatrices(element);
                endeixiStiffness = 2;
                //PrintUtilities.WriteToFile(Kt, @"C:\Users\turbo-x\Desktop\cohesive_check_MSOLVE_2\paradeigma_apo_arxika_swsta_shell_orthi_gia_check_tou_neou_class\orthi\CopyApoTaShellNewLoadCaseArgurhs\Kt_1.txt");
                IMatrix2D iGlobalStiffnessMatrix = new Matrix2D(Kt);
                //PrintUtilities.WriteToFile(Kt,@"C:\Users\turbo-x\Desktop\notes_elegxoi\MSOLVE_output\K_shell_arxiko.txt");
                return dofEnumerator.GetTransformedMatrix(iGlobalStiffnessMatrix);
            }
            else
            {
                //this.UpdatePartiallyPrecalculatedVariables_andforStiffnessMatrix(); // periexei this.updateck(); kai ta parakatw ean xreiasthei 
                //mporei na xwristhei kai na topothetithei edw to SPK_circumflex
                //Cons mporei na lamvanetai apo to uliko edw me MaterialsAtGP.ConstitutiveMatrix
                Kt = this.UpdateKmatrices(element);
                IMatrix2D iGlobalStiffnessMatrix = new Matrix2D(Kt);
                return dofEnumerator.GetTransformedMatrix(iGlobalStiffnessMatrix);
            }
        }

        public double[] CalculateForcesForLogging(Element element, double[] localDisplacements)
        {
            return CalculateForces(element, localDisplacements, new double[localDisplacements.Length]);
        }


    }




}









