using ISAAR.MSolve.FEM.Materials;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.Materials;
using ISAAR.MSolve.Materials.Interfaces; //using ISAAR.MSolve.PreProcessor.Interfaces;
using ISAAR.MSolve.MultiscaleAnalysis;
using ISAAR.MSolve.MultiscaleAnalysis.Interfaces;
using ISAAR.MSolve.MultiscaleAnalysis.SupportiveClasses;
using ISAAR.MSolve.MultiscaleAnalysisMerge;
using ISAAR.MSolve.Numerical.LinearAlgebra; //using ISAAR.MSolve.Matrices;
using ISAAR.MSolve.Solvers.Direct;
using ISAAR.MSolve.Solvers.Skyline;
using System;
using System.Collections.Generic;
using System.Linq;
using Xunit;

namespace ISAAR.MSolve.Tests.FEMpartB
{
    public class SeparateCodeCheckingClass4
    {              
        public static (double[], double[], double[], double[], IVector, IVector) Check05bStressIntegrationObje_v2_Integration()
        {
            //Origin: SeparateCodeCheckingClass.Check05bStressIntegration
            //modifications: tha xrhsimopoithei h nea microstructure me obje kapoia subdomainCalculations

            double E_disp = 3.5; /*Gpa*/ double ni_disp = 0.4; // stather Poisson
            ElasticMaterial3D material1 = new ElasticMaterial3D()
            { YoungModulus = E_disp, PoissonRatio = ni_disp, };
            double[,] DGtr = new double[3, 3] { { 1.10, 0, 0 }, { 0, 1, 0 }, { 0, 0, 1 } };
            double[] GLVec = Transform_DGtr_to_GLvec(DGtr);
            material1.UpdateMaterial(new StressStrainVectorContinuum3D(GLVec));
            //double[] stressesCheck1 = material1.Stresses;
            double[] stressesCheck1 = new double[6] {material1.Stresses[0], material1.Stresses[1], material1.Stresses[2],
                material1.Stresses[3],material1.Stresses[4],material1.Stresses[5] };
            DGtr = new double[3, 3] { { 1.20, 0, 0 }, { 0, 1, 0 }, { 0, 0, 1 } };
            GLVec = Transform_DGtr_to_GLvec(DGtr);
            material1.UpdateMaterial(new StressStrainVectorContinuum3D(GLVec));
            material1.SaveState();
            double[] stressesCheck2 = material1.Stresses.Data;

            //VectorExtensions.AssignTotalAffinityCount();
            IRVEbuilder_v2 homogeneousRveBuilder1 = new HomogeneousRVEBuilderNonLinear();
            //IRVEbuilder homogeneousRveBuilder1 = new HomogeneousRVEBuilderCheckEnaHexa();

            //IContinuumMaterial3DDefGrad_v2
            var microstructure3 = new MicrostructureDefGrad3D(homogeneousRveBuilder1, 
                model => (new SkylineSolver.Builder()).BuildSolver(model), false, 1);
            //IContinuumMaterial3DDefGrad microstructure3copyConsCheck = new Microstructure3copyConsCheckEna(homogeneousRveBuilder1);
            double[,] consCheck1 = new double[6, 6];
            for (int i1 = 0; i1 < 6; i1++) { for (int i2 = 0; i2 < 6; i2++) { consCheck1[i1, i2] = microstructure3.ConstitutiveMatrix[i1, i2]; } }

            microstructure3.UpdateMaterial(new double[9] { 1.10, 1, 1, 0, 0, 0, 0, 0, 0 });
            double[] stressesCheck3 = microstructure3.Stresses;
            microstructure3.SaveState();
            IVector uInitialFreeDOFs_state1 = microstructure3.uInitialFreeDOFDisplacementsPerSubdomain[1].Copy();

            microstructure3.UpdateMaterial(new double[9] { 1.20, 1, 1, 0, 0, 0, 0, 0, 0 });
            double[] stressesCheck4 = microstructure3.Stresses;
            IVector uInitialFreeDOFs_state2 = microstructure3.uInitialFreeDOFDisplacementsPerSubdomain[1].Copy();
            //PrintUtilities.WriteToFileVector(uInitialFreeDOFs_state1.CopyToArray(), @"C:\Users\turbo-x\Desktop\notes_elegxoi\MSOLVE_output_2\uInitialFreeDOFs_state1.txt");
            //PrintUtilities.WriteToFileVector(uInitialFreeDOFs_state2.CopyToArray(), @"C:\Users\turbo-x\Desktop\notes_elegxoi\MSOLVE_output_2\uInitialFreeDOFs_state2.txt");

            return (stressesCheck1, stressesCheck2, stressesCheck3, stressesCheck4, uInitialFreeDOFs_state1, uInitialFreeDOFs_state2);

        }

        public static (double[], double[], double[,], IVector, IVector) Check_Graphene_rve_Obje_v2_Integration()
        {
            //Origin: SeparateCodeCheckingClass4.Check05bStressIntegrationObje_v2_Integration parontos
            //modifications: tha xrhsimopoithei o GrapheneBuilder...35...v2 gia epilush enos paradeigmatos GrapheneReinforcedRVEBuilderCHECK
            //gia elegxo twn newn domwn
            //PROSOXH gia na elegxei kai h defterh iteration u_sunol_micro_2 prepei na valoume ston graphenebuilder Addgraphenesheet xwris to bondslip.

            double E_disp = 3.5; /*Gpa*/ double ni_disp = 0.4; // stather Poisson
            ElasticMaterial3D material1 = new ElasticMaterial3D()
            { YoungModulus = E_disp, PoissonRatio = ni_disp, };
            double[,] DGtr = new double[3, 3] { { 1.10, 0, 0 }, { 0, 1, 0 }, { 0, 0, 1 } };
            double[] GLVec = Transform_DGtr_to_GLvec(DGtr);
            material1.UpdateMaterial(new StressStrainVectorContinuum3D(GLVec));
            //double[] stressesCheck1 = material1.Stresses;
            double[] stressesCheck1 = new double[6] {material1.Stresses[0], material1.Stresses[1], material1.Stresses[2],
                material1.Stresses[3],material1.Stresses[4],material1.Stresses[5] };
            DGtr = new double[3, 3] { { 1.20, 0, 0 }, { 0, 1, 0 }, { 0, 0, 1 } };
            GLVec = Transform_DGtr_to_GLvec(DGtr);
            material1.UpdateMaterial(new StressStrainVectorContinuum3D(GLVec));
            material1.SaveState();
            double[] stressesCheck2 = material1.Stresses.Data;

            //VectorExtensions.AssignTotalAffinityCount();
            IRVEbuilder_v2 homogeneousRveBuilder1 = new RveGrShOne(1);
            //IRVEbuilder homogeneousRveBuilder1 = new HomogeneousRVEBuilderCheckEnaHexa();

            //IContinuumMaterial3DDefGrad_v2 
            var microstructure3 = new MicrostructureDefGrad3D(homogeneousRveBuilder1,
                model => (new SkylineSolver.Builder()).BuildSolver(model), false, 1);
            //IContinuumMaterial3DDefGrad microstructure3copyConsCheck = new Microstructure3copyConsCheckEna(homogeneousRveBuilder1);
            double[,] consCheck1 = new double[6, 6];
            for (int i1 = 0; i1 < 6; i1++) { for (int i2 = 0; i2 < 6; i2++) { consCheck1[i1, i2] = microstructure3.ConstitutiveMatrix[i1, i2]; } }

            microstructure3.UpdateMaterial(new double[9] { 1.05, 1, 1, 0, 0, 0, 0, 0, 0 });
            double[] stressesCheck3 = microstructure3.Stresses;
            microstructure3.SaveState();
            IVector uInitialFreeDOFs_state1 = microstructure3.uInitialFreeDOFDisplacementsPerSubdomain[1].Copy();

            microstructure3.UpdateMaterial(new double[9] { 1.10, 1, 1, 0, 0, 0, 0, 0, 0 });
            double[] stressesCheck4 = microstructure3.Stresses;
            IVector uInitialFreeDOFs_state2 = microstructure3.uInitialFreeDOFDisplacementsPerSubdomain[1].Copy();

            //PrintUtilities.WriteToFileVector(stressesCheck3, @"C:\Users\turbo-x\Desktop\notes_elegxoi\MSOLVE_output_2\stressesCheck3.txt");
            //PrintUtilities.WriteToFileVector(stressesCheck4, @"C:\Users\turbo-x\Desktop\notes_elegxoi\MSOLVE_output_2\stressesCheck4.txt");
            //PrintUtilities.WriteToFile(consCheck1, @"C:\Users\turbo-x\Desktop\notes_elegxoi\MSOLVE_output_2\consCheck1.txt");
            //PrintUtilities.WriteToFileVector(uInitialFreeDOFs_state1.CopyToArray(), @"C:\Users\turbo-x\Desktop\notes_elegxoi\MSOLVE_output_2\uInitialFreeDOFs_state1.txt");
            //PrintUtilities.WriteToFileVector(uInitialFreeDOFs_state2.CopyToArray(), @"C:\Users\turbo-x\Desktop\notes_elegxoi\MSOLVE_output_2\uInitialFreeDOFs_state2.txt");


            return (stressesCheck3, stressesCheck4, consCheck1, uInitialFreeDOFs_state1, uInitialFreeDOFs_state2);
        }

        #region transformation methods
        public static double[] Transform_DGtr_to_GLvec(double[,] DGtr)
        {
            double[,] GL = new double[3, 3];

            //
            for (int m = 0; m < 3; m++)
            {
                for (int n = 0; n < 3; n++)
                {
                    GL[m, n] = 0;
                    for (int p = 0; p < 3; p++)
                    {
                        GL[m, n] += DGtr[m, p] * DGtr[n, p];
                    }
                }
            }
            for (int m = 0; m < 3; m++)
            {
                GL[m, m] += -1;
            }
            for (int m = 0; m < 3; m++)
            {
                for (int n = 0; n < 3; n++)
                {
                    GL[m, n] = 0.5 * GL[m, n];
                }
            }

            double[] GLvec = new double[6];
            //
            for (int m = 0; m < 3; m++)
            {
                GLvec[m] = GL[m, m];
            }
            GLvec[3] = 2 * GL[0, 1];
            GLvec[4] = 2 * GL[1, 2];
            GLvec[5] = 2 * GL[2, 0];

            return GLvec;
        }
        #endregion
    }
}
