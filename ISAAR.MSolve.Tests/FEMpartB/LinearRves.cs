using ISAAR.MSolve.LinearAlgebra.Commons;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.Materials;
using ISAAR.MSolve.Materials.Interfaces;
using System;
using System.Collections.Generic;
using System.Text;
using Xunit;
using ISAAR.MSolve.MultiscaleAnalysis.Interfaces;
using ISAAR.MSolve.MultiscaleAnalysis;
using ISAAR.MSolve.Solvers.Direct;
using ISAAR.MSolve.FEM.Materials;

namespace ISAAR.MSolve.Tests.FEMpartB
{
    public static class LinearRves
    {
        [Fact]
        public static void CheckShellScaleTransitionsAndMicrostructure()
        {
            //Origin: Check05fStressIntegration
            //alllages: use of updated v2 classes


            double E_disp = 3.5; /*Gpa*/ double ni_disp = 0.4; // stather Poisson
            //var material1 = new ElasticMaterial2D(StressState2D.PlaneStress)
            //{ YoungModulus = E_disp, PoissonRatio = ni_disp, };
            //double[] GLVec = new double[3] { 0.01, 0, 0 };
            //material1.UpdateMaterial(new StressStrainVectorContinuum2D(GLVec));
            //double[] stressesCheck1 = new double[3] { material1.Stresses[0], material1.Stresses[1], material1.Stresses[2] };

            var Vec1 = Vector.CreateFromArray(new double[3] { 1, 0, 0 });
            var Vec2 = Vector.CreateFromArray(new double[3] { 0.5, 2, 0 });
            var strain = new double[3] { 0.01, 0, 0 };

            //var material2 = new ShellElasticMaterial2D() { YoungModulus = E_disp, PoissonRatio = ni_disp, TangentVectorV1 = new double[3] { Vec1[0], Vec1[1], Vec1[2] }, TangentVectorV2 = new double[3] { Vec2[0], Vec2[1], Vec2[2] } };
            //material2.UpdateMaterial(strain);
            //double[] stressesCheck2 = new double[3] { material2.Stresses[0], material2.Stresses[1], material2.Stresses[2] };

            var material3 = new ShellElasticMaterial2Dtransformationb() { YoungModulus = E_disp, PoissonRatio = ni_disp, TangentVectorV1 = new double[3] { Vec1[0], Vec1[1], Vec1[2] }, TangentVectorV2 = new double[3] { Vec2[0], Vec2[1], Vec2[2] } };
            var Matrix1 = Matrix.CreateZero(3, 3); for (int i1 = 0; i1 < 3; i1++) { for (int i2 = 0; i2 < 3; i2++) { Matrix1[i1, i2] = material3.ConstitutiveMatrix[i1, i2]; } }
            material3.UpdateMaterial(strain);
            double[] stressesCheck3 = new double[3] { material3.Stresses[0], material3.Stresses[1], material3.Stresses[2] };

            //VectorExtensions.AssignTotalAffinityCount();
            IdegenerateRVEbuilder_v2 homogeneousRveBuilder1 = new HomogeneousRVEBuilderLinearAndDegenerate();
            //IRVEbuilder homogeneousRveBuilder1 = new HomogeneousRVEBuilderCheckEnaHexa();

            //IContinuumMaterial2D microstructure3 = new Microstructure3DevelopMultipleSubdomainsUseBaseSmallStrains2D(homogeneousRveBuilder1);
            //IContinuumMaterial3DDefGrad microstructure3copyConsCheck = new Microstructure3copyConsCheckEna(homogeneousRveBuilder1);
            double[,] consCheck1 = new double[3, 3];
            for (int i1 = 0; i1 < 3; i1++) { for (int i2 = 0; i2 < 3; i2++) { consCheck1[i1, i2] = material3.ConstitutiveMatrix[i1, i2]; } }

            var material4 = new MicrostructureShell2D(homogeneousRveBuilder1, 
                model => (new SkylineSolver.Builder()).BuildSolver(model), false, 1)
            {
                TangentVectorV1 = new double[3] { Vec1[0], Vec1[1], Vec1[2] },
                TangentVectorV2 = new double[3] { Vec2[0], Vec2[1], Vec2[2] }
            };
            var Matrix2 = Matrix.CreateZero(3, 3); for (int i1 = 0; i1 < 3; i1++) { for (int i2 = 0; i2 < 3; i2++) { Matrix2[i1, i2] = material4.ConstitutiveMatrix[i1, i2]; } }
            material4.UpdateMaterial(strain);
            double[] stressesCheck4 = new double[3] { material4.Stresses[0], material4.Stresses[1], material4.Stresses[2] };


            //-------------Check 2 steps savestate etc---------------
            material4.SaveState();
            material4.UpdateMaterial(new double[3] { 2 * strain[0], 2 * strain[1], 2 * strain[2] });
            double[] stressesCheck5 = new double[3] { material4.Stresses[0], material4.Stresses[1], material4.Stresses[2] };

            Assert.True(NRNLAnalyzerDevelopTest_v2.AreDisplacementsSame_v2(stressesCheck3, stressesCheck4));
            Assert.True(NRNLAnalyzerDevelopTest_v2.AreDisplacementsSame_v2(new double[3] { 2 * stressesCheck3[0], 2 * stressesCheck3[1], 2 * stressesCheck3[2] },
                                                                            stressesCheck5));
            Assert.True(BondSlipTest.AreDisplacementsSame_v2(Matrix1.CopyToArray2D(), consCheck1));
            Assert.True(AreDisplacementsSame_v2(Matrix1.CopyToArray2D(), material4.ConstitutiveMatrix));
        }

        [Fact]
        public static void Check2DscaleTransitionsAndMicrostructure()
        {        
            //Check05cStressIntegration()
            double E_disp = 3.5; /*Gpa*/ double ni_disp = 0.4; // stathera Poisson
            var material1 = new ElasticMaterial2D_v2 (StressState2D.PlaneStress)
            { YoungModulus = E_disp, PoissonRatio = ni_disp, };
            double[] GLVec = new double[3] { 0.01, 0, 0 };
            material1.UpdateMaterial(GLVec);
            double[] stressesCheck1 = new double[3] { material1.Stresses[0], material1.Stresses[1], material1.Stresses[2] };
            //material1.SaveState();
            GLVec = new double[3] { 0.02, 0, 0 };
            material1.UpdateMaterial(GLVec);
            double[] stressesCheck2 = new double[3] { material1.Stresses[0], material1.Stresses[1], material1.Stresses[2] };

            //VectorExtensions.AssignTotalAffinityCount();
            IdegenerateRVEbuilder_v2 homogeneousRveBuilder1 = new HomogeneousRVEBuilderLinearAndDegenerate();
            //IRVEbuilder homogeneousRveBuilder1 = new HomogeneousRVEBuilderCheckEnaHexa();

            IContinuumMaterial2D_v2 microstructure3 = new Microstructure2DplaneStress(homogeneousRveBuilder1, 
                model => (new SkylineSolver.Builder()).BuildSolver(model), false, 1);
            //IContinuumMaterial3DDefGrad microstructure3copyConsCheck = new Microstructure3copyConsCheckEna(homogeneousRveBuilder1);
            double[,] consCheck1 = new double[3, 3];
            for (int i1 = 0; i1 < 3; i1++) { for (int i2 = 0; i2 < 3; i2++) { consCheck1[i1, i2] = microstructure3.ConstitutiveMatrix[i1, i2]; } }

            microstructure3.UpdateMaterial(new double[3] { 0.010, 0, 0 });
            double[] stressesCheck3 = new double[3] { microstructure3.Stresses[0], microstructure3.Stresses[1], microstructure3.Stresses[2] };
            microstructure3.SaveState();
            microstructure3.UpdateMaterial(new double[3] { 0.020, 0, 0 });
            double[] stressesCheck4 = new double[3] { microstructure3.Stresses[0], microstructure3.Stresses[1], microstructure3.Stresses[2] };

            microstructure3.SaveState();
            microstructure3.UpdateMaterial(new double[3] { 0.030, 0, 0 });
            double[] stressesCheck5 = new double[3] { microstructure3.Stresses[0], microstructure3.Stresses[1], microstructure3.Stresses[2] };
            var Matrix1 = Matrix.CreateZero(3, 3); for (int i1 = 0; i1 < 3; i1++) { for (int i2 = 0; i2 < 3; i2++) { Matrix1[i1, i2] = microstructure3.ConstitutiveMatrix[i1, i2]; } }

            Assert.True(NRNLAnalyzerDevelopTest_v2.AreDisplacementsSame_v2(stressesCheck1, stressesCheck3));
            Assert.True(NRNLAnalyzerDevelopTest_v2.AreDisplacementsSame_v2(stressesCheck2, stressesCheck4));
            Assert.True(NRNLAnalyzerDevelopTest_v2.AreDisplacementsSame_v2(new double[3] { 3 * stressesCheck1[0], 3 * stressesCheck1[1], 3 * stressesCheck1[2] },
                                                                            stressesCheck5));
            Assert.True(AreDisplacementsSame_v2( consCheck1, material1.ConstitutiveMatrix));
            Assert.True(AreDisplacementsSame_v2(Matrix1.CopyToArray2D(), material1.ConstitutiveMatrix));
        }

        [Fact]
        public static void Check3DscaleTransitionsAndMicrostructure()
        {
            //Check05c2_3D_StressIntegration
            double E_disp = 3.5; /*Gpa*/ double ni_disp = 0.4; // stather Poisson
            var material1 = new ElasticMaterial3D_v2()
            { YoungModulus = E_disp, PoissonRatio = ni_disp, };
            double[] GLVec = new double[6] { 0.01, 0, 0, 0, 0, 0 };
            material1.UpdateMaterial(GLVec);
            double[] stressesCheck1 = new double[6] { material1.Stresses[0], material1.Stresses[1], material1.Stresses[2], material1.Stresses[3], material1.Stresses[4], material1.Stresses[5] };
            //material1.SaveState();
            GLVec = new double[6] { 0, 0, 0, 0, 0.02, 0 };
            material1.UpdateMaterial(GLVec);
            double[] stressesCheck2 = new double[6] { material1.Stresses[0], material1.Stresses[1], material1.Stresses[2], material1.Stresses[3], material1.Stresses[4], material1.Stresses[5] };

            //VectorExtensions.AssignTotalAffinityCount();
            IRVEbuilder_v2 homogeneousRveBuilder1 = new HomogeneousRVEBuilderLinear();
            //IRVEbuilder homogeneousRveBuilder1 = new HomogeneousRVEBuilderCheckEnaHexa();

            IContinuumMaterial3D_v2 microstructure3 = new Microstructure3D(homogeneousRveBuilder1, 
                model => (new SkylineSolver.Builder()).BuildSolver(model), false, 1);
            //IContinuumMaterial3DDefGrad microstructure3copyConsCheck = new Microstructure3copyConsCheckEna(homogeneousRveBuilder1);
            double[,] consCheck1 = new double[6, 6];
            for (int i1 = 0; i1 < 6; i1++) { for (int i2 = 0; i2 < 6; i2++) { consCheck1[i1, i2] = microstructure3.ConstitutiveMatrix[i1, i2]; } }

            microstructure3.UpdateMaterial(new double[6] { 0.010, 0, 0, 0, 0, 0 });
            double[] stressesCheck3 = new double[6] { microstructure3.Stresses[0], microstructure3.Stresses[1], microstructure3.Stresses[2], microstructure3.Stresses[3], microstructure3.Stresses[4], microstructure3.Stresses[5] };
            microstructure3.SaveState();
            microstructure3.UpdateMaterial(new double[6] { 0, 0, 0, 0, 0.020, 0 });
            double[] stressesCheck4 = new double[6] { microstructure3.Stresses[0], microstructure3.Stresses[1], microstructure3.Stresses[2], microstructure3.Stresses[3], microstructure3.Stresses[4], microstructure3.Stresses[5] };

            microstructure3.SaveState();
            microstructure3.UpdateMaterial(new double[6] { 0.030, 0, 0, 0, 0, 0 });
            double[] stressesCheck5 = new double[6] { microstructure3.Stresses[0], microstructure3.Stresses[1], microstructure3.Stresses[2], microstructure3.Stresses[3], microstructure3.Stresses[4], microstructure3.Stresses[5] };
            var Matrix1 = Matrix.CreateZero(3, 3); for (int i1 = 0; i1 < 3; i1++) { for (int i2 = 0; i2 < 3; i2++) { Matrix1[i1, i2] = microstructure3.ConstitutiveMatrix[i1, i2]; } }

            Assert.True(NRNLAnalyzerDevelopTest_v2.AreDisplacementsSame_v2(stressesCheck1, stressesCheck3));
            Assert.True(NRNLAnalyzerDevelopTest_v2.AreDisplacementsSame_v2(stressesCheck2, stressesCheck4));
            Assert.True(NRNLAnalyzerDevelopTest_v2.AreDisplacementsSame_v2(new double[6] { 3 * stressesCheck1[0], 3 * stressesCheck1[1], 3 * stressesCheck1[2], 3 * stressesCheck1[3], 3 * stressesCheck1[4], 3 * stressesCheck1[5] },
                                                                            stressesCheck5));
            Assert.True(AreDisplacementsSame_v2(consCheck1, material1.ConstitutiveMatrix));
            Assert.True(AreDisplacementsSame_v2(Matrix1.CopyToArray2D(), material1.ConstitutiveMatrix));
        }

        public static bool AreDisplacementsSame_v2(double[,] expectedValues,
            IMatrixView computedValues)
        {
            var comparer = new ValueComparer(1E-14);
            for (int i1 = 0; i1 < expectedValues.GetLength(0); i1++)
            {
                for (int i2 = 0; i2 < expectedValues.GetLength(1); i2++)
                {
                    if (!comparer.AreEqual(expectedValues[i1, i2], computedValues[i1, i2]))
                    {
                        return false;
                    }
                }
            }
            return true;
        }

    }
}
