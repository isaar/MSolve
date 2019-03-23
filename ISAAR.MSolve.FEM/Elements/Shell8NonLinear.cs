using ISAAR.MSolve.Discretization;
using ISAAR.MSolve.Discretization.Integration.Quadratures;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.FEM.Elements.SupportiveClasses;
using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.FEM.Interfaces;
using ISAAR.MSolve.FEM.Interpolation;
using ISAAR.MSolve.FEM.Interpolation.Jacobians;
using ISAAR.MSolve.Materials.Interfaces;
using ISAAR.MSolve.Numerical.LinearAlgebra;
using ISAAR.MSolve.Numerical.LinearAlgebra.Interfaces;
using System;
using System.Collections.Generic;

//TODO: move stuff to Shell8DirectionVectorUtilities
namespace ISAAR.MSolve.FEM.Elements
{
    public class Shell8NonLinear : IStructuralFiniteElement
    {
        protected readonly static DOFType[] nodalDOFTypes = new DOFType[] { DOFType.X, DOFType.Y, DOFType.Z, DOFType.RotX, DOFType.RotY };
        protected readonly static DOFType[][] dofTypes = new DOFType[][] { nodalDOFTypes, nodalDOFTypes, nodalDOFTypes,
            nodalDOFTypes, nodalDOFTypes, nodalDOFTypes, nodalDOFTypes, nodalDOFTypes };
        protected readonly IShellMaterial[] materialsAtGaussPoints; 
        protected IElementDOFEnumerator dofEnumerator = new GenericDOFEnumerator();
        
        public double[][] oVn_i { get; set; }
        public double[] tk { get; set; } 
        private int nGaussPoints;

        private double[] integrationCoefficient;
        private double[][] tU;   //dimensions 8 arrays of six elements
        private double[][] tUvec;//dimensions 8 arrays of six elements
        
        // updating element configurations (NESSESARY for UpdateCoordinateData method)
        // auxiliary fields for calculating element iterative rotations
        private double[] ak_total = new double[8];
        private double[] bk_total = new double[8];

        private double[][] GLvec; // TODO possibly gl_vec_last_converged can be saved too if strains aren't handled in tthe current way by shell material classes
        
        public Shell8NonLinear(IShellMaterial material, IQuadrature3D quadratureForStiffness)
        {
            this.Interpolation = InterpolationShell8.UniqueInstance;
            this.QuadratureForStiffness = quadratureForStiffness;
            this.nGaussPoints = quadratureForStiffness.IntegrationPoints.Count;
            materialsAtGaussPoints = new IShellMaterial[nGaussPoints];
            for (int i = 0; i < nGaussPoints; i++) materialsAtGaussPoints[i] = material.Clone();
        }

        public InterpolationShell8 Interpolation { get; }
        public IQuadrature3D QuadratureForStiffness { get; }
        
        private Matrix2D[] GetBL11a(Matrix2D[] J_0inv)
        {
            Matrix2D[] BL11a;
            BL11a = new Matrix2D[nGaussPoints];
            for (int j = 0; j < nGaussPoints; j++)
            { BL11a[j] = new Matrix2D(6, 9); }
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

                //calculating [4,4]:[4,6] and [5,7]:[5,9]
                for (int k = 0; k < 2; k++)
                {
                    for (int l = 0; l < 3; l++)
                    { BL11a[j][3 + k, 3 + 3 * k + l] = J_0inv[j][k, l]; }
                }

                //calculating  [4,1]:[4,3] and [5,4]:[5,6]
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

        private Matrix2D[] GetBL12(Matrix2D[] J_0inv)
        {
            Matrix2D[] BL12;
            BL12 = new Matrix2D[nGaussPoints];
            for (int j = 0; j < nGaussPoints; j++)
            { BL12[j] = new Matrix2D( 9, 9); }
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

        private Matrix2D[] GetBL13(IReadOnlyList<Matrix2D> shapeFunctionDerivatives,double [][] tUvec , Matrix2D[] J_0a)
        {
            Matrix2D[] BL13;
            BL13 = new Matrix2D[nGaussPoints];
            for (int j = 0; j < nGaussPoints; j++)
            { BL13[j] = new Matrix2D(9, 40); }
            for (int j = 0; j < nGaussPoints; j++)
            {
               
                //columns 1:3
                for (int m = 0; m < 8; m++)
                {
                    for (int k = 0; k < 3; k++)
                    {
                        for (int l = 0; l < 2; l++)
                        {
                            BL13[j][3 * k + l, 5 * m + k] = shapeFunctionDerivatives[j][m, l];
                        }
                    }
                }

                //columns 4:5
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

        private Matrix2D[] GetBNL1(Matrix2D[] J_0inv)
        {
            Matrix2D[] BNL1;
            BNL1 = new Matrix2D[nGaussPoints];
            for (int j = 0; j < nGaussPoints; j++)
            { BNL1[j] = new Matrix2D(9, 9); }
            for (int j = 0; j < nGaussPoints; j++)
            {
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
            double[] E;
            double[] ni;


            IReadOnlyList<Vector> shapeFunctions = Interpolation.EvaluateFunctionsAtGaussPoints(QuadratureForStiffness);
            IReadOnlyList<Matrix2D> shapeFunctionsDerivatives = Interpolation.EvaluateNaturalGradientsAtGaussPoints(QuadratureForStiffness);            
            (Matrix2D[] ll1, Matrix2D[] J_0a) = JacobianShell8Calculations.Getll1AndJ_0a(QuadratureForStiffness, tk, shapeFunctions, shapeFunctionsDerivatives);
            //TODO J_0, J_0b etc. can be cached for not large scale models

            tx_i = new double[8][];
            (tU, tUvec) = Shell8DirectionVectorUtilities.GetInitialDirectionVectorValues(oVn_i);
            double[][] oV1_i = new double[8][]; //tangent vector ''1'' initial configuration
            for (int j = 0; j < 8; j++)
            {
                tx_i[j] = new double[] { element.INodes[j].X, element.INodes[j].Y, element.INodes[j].Z, };
                oV1_i[j] = new double[3];

                oV1_i[j][0] = tUvec[j][0];
                oV1_i[j][1] = tUvec[j][1];
                oV1_i[j][2] = tUvec[j][2];

            }

            (Matrix2D[] J_0inv, double[] detJ_0) =
                JacobianShell8Calculations.GetJ_0invAndDetJ_0(J_0a, element.INodes, oVn_i, nGaussPoints);

            for (int j = 0; j < nGaussPoints; j++)
            {
                double[] V3 = new double[3];
                double V3_norm;
                double[] V1 = new double[3];
                double V1_norm;
                
                
                for (int k = 0; k < 3; k++)
                { V3[k] = 0; V1[k] = 0; /* V2[k] = 0; */ }

                for (int k = 0; k < 8; k++)
                {
                    for (int l = 0; l < 3; l++)
                    {
                        V3[l] += shapeFunctions[j][k] * oVn_i[k][l];
                        V1[l] += shapeFunctions[j][k] * oV1_i[k][l];
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

            }            

            integrationCoefficient = new double[nGaussPoints];
            for (int j = 0; j < nGaussPoints; j++)
            { integrationCoefficient[j] += QuadratureForStiffness.IntegrationPoints[j].Weight * detJ_0[j]; }

        }

       

        private double [][,] CalculateCk(Matrix2D[] J_0a,double [][] tU )
        {
            double[][,] ck;
            ck = new double[nGaussPoints][,];
            for (int j = 0; j < nGaussPoints; j++)
            {
                ck[j] = new double[8, 9];
            }
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

        private void CalculateStrains(IElement element, double[][] tx_i)
        {
            IReadOnlyList<Vector> shapeFunctions = Interpolation.EvaluateFunctionsAtGaussPoints(QuadratureForStiffness);
            IReadOnlyList<Matrix2D> shapeFunctionsDerivatives = Interpolation.EvaluateNaturalGradientsAtGaussPoints(QuadratureForStiffness);
            (Matrix2D[] ll1, Matrix2D[] J_0a) = JacobianShell8Calculations.Getll1AndJ_0a(QuadratureForStiffness, tk, shapeFunctions, shapeFunctionsDerivatives);

            (Matrix2D[] J_0inv, double[] detJ_0) =
                JacobianShell8Calculations.GetJ_0invAndDetJ_0(J_0a, element.INodes, oVn_i, nGaussPoints);

            Matrix2D[] J_1 = JacobianShell8Calculations.Get_J_1(nGaussPoints, tx_i, tU, J_0a);
            Matrix2D[] DefGradTr;
            Matrix2D[] GL;
            
            DefGradTr = new Matrix2D[nGaussPoints];
            GL = new Matrix2D[nGaussPoints];
            for (int j = 0; j < nGaussPoints; j++)
            {
                
                DefGradTr[j] = new Matrix2D(3, 3);
                GL[j] = new Matrix2D(3, 3);
            }

            GLvec = new double[nGaussPoints][]; 

            for (int j = 0; j < nGaussPoints; j++)
            {
                DefGradTr[j] = J_0inv[j] * J_1[j];
                
                GL[j] = DefGradTr[j] * DefGradTr[j].Transpose();
                
                for (int k = 0; k < 3; k++)
                {
                    GL[j][k, k] = GL[j][k, k] - 1;
                }
                GL[j].Scale(0.5);                

                GLvec[j] = new double[6];
                for (int k = 0; k < 3; k++)
                { GLvec[j][k] = GL[j][k, k]; }
                GLvec[j][3] = 2 * GL[j][0, 1];
                GLvec[j][4] = 2 * GL[j][1, 2];
                GLvec[j][5] = 2 * GL[j][2, 0];

            }
        }               

        private double [] UpdateForces(IElement element)
        {

            IReadOnlyList<Vector> shapeFunctions = Interpolation.EvaluateFunctionsAtGaussPoints(QuadratureForStiffness);
            IReadOnlyList<Matrix2D> shapeFunctionsDerivatives = Interpolation.EvaluateNaturalGradientsAtGaussPoints(QuadratureForStiffness);
            (Matrix2D[] ll1, Matrix2D[] J_0a) = JacobianShell8Calculations.Getll1AndJ_0a(QuadratureForStiffness, tk, shapeFunctions, shapeFunctionsDerivatives);
            
            //BL11a etc. are not cached currently(declare here)
            (Matrix2D[] J_0inv, double[] detJ_0) =
                JacobianShell8Calculations.GetJ_0invAndDetJ_0(J_0a, element.INodes, oVn_i, nGaussPoints);
            Matrix2D[] BL11a;
            BL11a = GetBL11a(J_0inv);
            Matrix2D[] BL12;
            BL12 = GetBL12(J_0inv);
            Matrix2D[] BL13;
            BL13 = GetBL13(shapeFunctionsDerivatives, tUvec, J_0a);

            Matrix2D[] BL;
            BL = new Matrix2D[nGaussPoints];
            Matrix2D[] BL01plus1_2;
            BL01plus1_2 = new Matrix2D[nGaussPoints];
            for (int j = 0; j < nGaussPoints; j++)
            {
                BL[j] = new Matrix2D(6, 40);
                BL01plus1_2[j] = new Matrix2D(6, 9);
            }


            double[][] Fxk= new double [nGaussPoints+1] [];
            for (int j = 0; j<nGaussPoints + 1; j++)
            {
                Fxk[j] = new double[40];
            }

            Matrix2D ll2;
            ll2 = new Matrix2D(24, 3);                 
            for (int j = 0; j < 8; j++)
            {
                for (int k = 0; k < 3; k++)
                {
                    ll2[3 * j + 0, k] = tU[j][k];
                    ll2[3 * j + 1, k] = tU[j][3 + k];
                    ll2[3 * j + 2, k] = oVn_i[j][k];
                }
            }

            Matrix2D[] l_circumflex;
            l_circumflex = new Matrix2D[nGaussPoints];
            for (int j = 0; j < nGaussPoints; j++)
            { l_circumflex[j] = new Matrix2D(3, 3); }
            for (int j = 0; j < nGaussPoints; j++)
            {
                l_circumflex[j] = ll1[j] * ll2;                
            }

            Matrix2D[] BL11b;
            BL11b = new Matrix2D[nGaussPoints];
            for (int j = 0; j < nGaussPoints; j++)
            { BL11b[j] = new Matrix2D(9, 9); }
            for (int j = 0; j < nGaussPoints; j++)
            {
                
                for (int k = 0; k < 3; k++)
                {
                    for (int l = 0; l < 3; l++)
                    {
                        for (int m = 0; m < 3; m++)
                        { BL11b[j][3 * k + l, 3 * k + m] = l_circumflex[j][l, m]; } //don;t change that
                    }
                }
            }

            Matrix2D[] BL11;
            BL11 = new Matrix2D[nGaussPoints];            
            for (int j = 0; j < nGaussPoints; j++)
            {
                BL11[j] = BL11a[j] * BL11b[j];           
            }

            
            Matrix2D[] BL1_2;
            BL1_2 = new Matrix2D[nGaussPoints];
            for (int j = 0; j < nGaussPoints; j++)
            {
                BL1_2[j] = BL11[j] * BL12[j];
                
                for (int k = 0; k < 6; k++)
                {
                    for (int l = 0; l < 9; l++)
                    {
                        BL01plus1_2[j][k, l] = BL1_2[j][k, l] + BL11a[j][k, l];
                    }
                }

                BL[j] = BL01plus1_2[j] * BL13[j];                 
            }

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

        private Matrix2D UpdateKmatrices(IElement element)
        {            
            double[][] kck;// 1 per node and (one per Gauss point+1 for the addition) -->[GP][8 dofs_per_node*nodes]
            Matrix2D[] KNL;
            Matrix2D[] KL;
            double[][] BL01plus1_2tSPKvec;


            IReadOnlyList<Vector> shapeFunctions = Interpolation.EvaluateFunctionsAtGaussPoints(QuadratureForStiffness);
            IReadOnlyList<Matrix2D> shapeFunctionsDerivatives = Interpolation.EvaluateNaturalGradientsAtGaussPoints(QuadratureForStiffness);
            (Matrix2D[] ll1, Matrix2D[] J_0a) = JacobianShell8Calculations.Getll1AndJ_0a(QuadratureForStiffness, tk, shapeFunctions, shapeFunctionsDerivatives);
            (Matrix2D[] J_0inv, double[] detJ_0) =
                JacobianShell8Calculations.GetJ_0invAndDetJ_0(J_0a, element.INodes, oVn_i, nGaussPoints);
            Matrix2D[] BNL1;
            BNL1 = GetBNL1(J_0inv);
            Matrix2D[] BL13;
            BL13 = GetBL13(shapeFunctionsDerivatives, tUvec, J_0a); //TODO: maybe cached from calcForces for problems of normal memory requirements
            double[][,] ck;
            ck = CalculateCk(J_0a, tU);

            //
            Matrix2D[] BL11a; //TODO: maybe cached from calcForces for problems of normal memory
            BL11a = GetBL11a(J_0inv);
            Matrix2D[] BL12;
            BL12 = GetBL12(J_0inv);
            //
            Matrix2D[] BL; //TODO: maybe cached from calcForces for problems of normal memory
            BL = new Matrix2D[nGaussPoints];
            Matrix2D[] BL01plus1_2;
            BL01plus1_2 = new Matrix2D[nGaussPoints];
            for (int j = 0; j < nGaussPoints; j++)
            {
                BL[j] = new Matrix2D(6, 40);
                BL01plus1_2[j] = new Matrix2D(6, 9);
            }
            //
            var ll2 = new Matrix2D(24, 3);
            for (int j = 0; j < 8; j++)
            {
                for (int k = 0; k < 3; k++)
                {
                    ll2[3 * j + 0, k] = tU[j][k];
                    ll2[3 * j + 1, k] = tU[j][3 + k];
                    ll2[3 * j + 2, k] = oVn_i[j][k];
                }
            }

            Matrix2D[] l_circumflex;
            l_circumflex = new Matrix2D[nGaussPoints];
            for (int j = 0; j < nGaussPoints; j++)
            { l_circumflex[j] = new Matrix2D(3, 3); }
            for (int j = 0; j < nGaussPoints; j++)
            {
                l_circumflex[j] = ll1[j] * ll2;                           
            }

            Matrix2D[] BL11b;
            BL11b = new Matrix2D[nGaussPoints];
            for (int j = 0; j < nGaussPoints; j++)
            { BL11b[j] = new Matrix2D(9, 9); }
            for (int j = 0; j < nGaussPoints; j++)
            {
                for (int k = 0; k < 3; k++)
                {
                    for (int l = 0; l < 3; l++)
                    {
                        for (int m = 0; m < 3; m++)
                        { BL11b[j][3 * k + l, 3 * k + m] = l_circumflex[j][l, m]; }
                    }
                }
            }

            Matrix2D[] BL11;
            BL11 = new Matrix2D[nGaussPoints];
            for (int j = 0; j < nGaussPoints; j++)
            { BL11[j] = new Matrix2D(6, 9); }
            for (int j = 0; j < nGaussPoints; j++)
            {
                BL11[j] = BL11a[j] * BL11b[j];                
            }


            Matrix2D[] BL1_2;
            BL1_2 = new Matrix2D[nGaussPoints];
            for (int j = 0; j < nGaussPoints; j++)
            {
                BL1_2[j] = new Matrix2D(6, 9);
            }
            for (int j = 0; j < nGaussPoints; j++)
            {
                BL1_2[j] = BL11[j] * BL12[j];
                
                for (int k = 0; k < 6; k++)
                {
                    for (int l = 0; l < 9; l++)
                    {
                        BL01plus1_2[j][k, l] = BL1_2[j][k, l] + BL11a[j][k, l];
                    }
                }

                BL[j] = BL01plus1_2[j] * BL13[j];               
            }

            KL= new Matrix2D[nGaussPoints + 1];
            KNL = new Matrix2D[nGaussPoints + 1];
            kck = new double[nGaussPoints + 1][];
            BL01plus1_2tSPKvec = new double[nGaussPoints][];
            for (int j = 0; j < nGaussPoints + 1; j++)
            {
                
                KL[j] = new Matrix2D(40, 40);
                KNL[j] = new Matrix2D(40, 40);
            }
            
           for (int j = 0; j < nGaussPoints; j++)
            {
                var SPK_circumflex = new Matrix2D(9, 9);                
                for (int k = 0; k < 3; k++) 
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

                var BNL = BNL1[j] * BL13[j];

                BL01plus1_2tSPKvec[j] = new double[9];
                for (int k = 0; k < 9; k++)
                {
                    for (int m = 0; m < 6; m++)
                    {
                        BL01plus1_2tSPKvec[j][k] += BL01plus1_2[j][m, k] * materialsAtGaussPoints[j].Stresses[m];
                    }
                }

                kck[j] = new double[8];
                for (int k = 0; k < 8; k++)
                {
                    for (int m = 0; m < 9; m++)
                    {
                        kck[j][k] += ck[j][k, m] * BL01plus1_2tSPKvec[j][m];
                    }
                }

                
                var ConsBL = new Matrix2D(6, 40);
                for (int k = 0; k < 6; k++)
                {
                    for (int l = 0; l < 40; l++)
                    {                        
                        for (int m = 0; m < 6; m++)
                        {
                            ConsBL[k, l] += materialsAtGaussPoints[j].ConstitutiveMatrix[k, m] * BL[j][m, l];
                        }
                    }
                }

                var S_BNL = SPK_circumflex * BNL;                
                KNL[j] = BNL.Transpose() * S_BNL;                
                KL[j] = BL[j].Transpose() * ConsBL;
                
            }

            var Kt = new Matrix2D(40, 40);
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

        private void UpdateCoordinateData(double[] localdisplacements, IList<INode> elementNodes, out double[][] tx_i)
        {
            double[][] ox_i = new double[8][]; // this should be allocated in the constructor
            for (int j = 0; j < 8; j++)
            {
                ox_i[j] = new double[] { elementNodes[j].X, elementNodes[j].Y, elementNodes[j].Z, };
            }

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
                //update tU andUvec 
                tU[k][0] = localdisplacements[5 * k + 0];
                tU[k][1] = localdisplacements[5 * k + 1];
                tU[k][2] = localdisplacements[5 * k + 2];
                ak = localdisplacements[5 * k + 3] - ak_total[k];
                ak_total[k] = localdisplacements[5 * k + 3];
                bk = localdisplacements[5 * k + 4] - bk_total[k];
                bk_total[k] = localdisplacements[5 * k + 4];
                Shell8DirectionVectorUtilities.RotateNodalDirectionVectors(ak, bk, k,tU,tUvec);                
            }            
        }
       
        // region IstructuralFiniteElement

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

        //IstructuralElement concerning material
        public void ClearMaterialState()
        {
            foreach (IShellMaterial m in materialsAtGaussPoints) m.ClearState();
        }
        public void SaveMaterialState()
        {
            foreach (IShellMaterial m in materialsAtGaussPoints) m.SaveState();
        }

        public void ClearMaterialStresses()
        {
            foreach (IShellMaterial m in materialsAtGaussPoints) m.ClearStresses();
        }

        public bool MaterialModified
        {
            get
            {
                foreach (IShellMaterial material in materialsAtGaussPoints)
                    if (material.Modified) return true;
                return false;
            }
        }

        public void ResetMaterialModified()
        {
            foreach (IShellMaterial material in materialsAtGaussPoints) material.ResetModified();
        }

        public Tuple<double[], double[]> CalculateStresses(Element element, double[] localTotalDisplacements, double[] localdDisplacements)
        {
            this.UpdateCoordinateData(localTotalDisplacements, element.INodes, out double[][] tx_i);
            this.CalculateStrains(element, tx_i);
            for (int npoint = 0; npoint < materialsAtGaussPoints.Length; npoint++)
            {
                materialsAtGaussPoints[npoint].UpdateMaterial(GLvec[npoint]);
            }
           
            return new Tuple<double[], double[]>(new double[123], materialsAtGaussPoints[materialsAtGaussPoints.Length - 1].Stresses);
        }

        //Istructural: dynamic
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

        public double[] CalculateForces(Element element, double[] localTotalDisplacements, double[] localdDisplacements)
        {
            return UpdateForces(element);
        }

        private int endeixiStiffness = 1;
        public virtual IMatrix2D StiffnessMatrix(IElement element)
        {
            var Kt = new Matrix2D(40, 40);
            if (endeixiStiffness == 1)
            {
                this.CalculateInitialConfigurationData(element, out double[][] tx_i);
                this.CalculateStrains(element, tx_i);

                for (int npoint = 0; npoint < materialsAtGaussPoints.Length; npoint++)
                {
                    materialsAtGaussPoints[npoint].UpdateMaterial(GLvec[npoint]);
                }

                this.UpdateForces(element); 


                Kt = this.UpdateKmatrices(element);
                endeixiStiffness = 2;
                return dofEnumerator.GetTransformedMatrix(Kt);
            }
            else
            {
                Kt = this.UpdateKmatrices(element);
                return dofEnumerator.GetTransformedMatrix(Kt);
            }
        }

        public double[] CalculateForcesForLogging(Element element, double[] localDisplacements)
        {
            return CalculateForces(element, localDisplacements, new double[localDisplacements.Length]);
        }


    }




}









