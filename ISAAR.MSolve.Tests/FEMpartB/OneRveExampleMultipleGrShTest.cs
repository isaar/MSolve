using ISAAR.MSolve.LinearAlgebra.Commons;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.Materials;
using ISAAR.MSolve.Materials.Interfaces;
using ISAAR.MSolve.MultiscaleAnalysis.SupportiveClasses;
using ISAAR.MSolve.Tests.FEMpartB.SeparationBenchmarks2;
using System;
using System.Collections.Generic;
using System.Text;
using Xunit;

namespace ISAAR.MSolve.Tests.FEMpartB
{
    public static class OneRveExampleMultipleGrShTest
    {
        [Fact]
        public static void CheckOneRveSerial()
        {
            (double[] stressesCheck3, double[] stressesCheck4, double[,] consCheck1, IVector uInitialFreeDOFs_state1, IVector uInitialFreeDOFs_state2) = OneRveExample.Check_Graphene_rve_serial();

            string results_file1 = "..\\..\\..\\InputFiles\\MicroStructureAndIntegrationOneRveMultipleGrSh\\uInitialFreeDOFs_state1.txt";
            string results_file2 = "..\\..\\..\\InputFiles\\MicroStructureAndIntegrationOneRveMultipleGrSh\\uInitialFreeDOFs_state2.txt";
            //string results_file3 = "..\\..\\..\\InputFiles\\MicroStructureAndIntegrationTestsB\\consCheck1.txt";
            string results_file4 = "..\\..\\..\\InputFiles\\MicroStructureAndIntegrationOneRveMultipleGrSh\\stressesCheck3.txt";
            string results_file5 = "..\\..\\..\\InputFiles\\MicroStructureAndIntegrationOneRveMultipleGrSh\\stressesCheck4.txt";
            double[] displacements1sIncrement = PrintUtilities.ReadVector(results_file1);
            double[] displacements2ndIncrement = PrintUtilities.ReadVector(results_file2);

            double[,] consCheck1Expected = new double[6, 6]
            {{8.247794441602281,5.040483382644040,5.045179342838760,-0.034573545680066,0.012873618640199,0.067413461733790},
            {5.040483382644040,7.758675250745090,5.083447516662590,-0.017660393516958,0.086264761000810,-0.001886483315119},
            {5.045179342838760,5.083447516662600,7.889514025249530,0.014993568822868,0.174547712576532,0.013639601528685},
            {-0.034573545680067,-0.017660393516956,0.014993568822868,1.404689076704550,0.023343385610862,0.099337624448147},
            {0.012873618640199,0.086264761000810,0.174547712576533,0.023343385610861,1.347276707954930,-0.002212957880199},
            {0.067413461733791,-0.001886483315119,0.013639601528686,0.099337624448147,-0.002212957880199,1.454060010268960} };

            double[] stressesCheck3Expected = PrintUtilities.ReadVector(results_file4);
            double[] stressesCheck4Expected = PrintUtilities.ReadVector(results_file5);

            Assert.True(NRNLAnalyzerDevelopTest_v2.AreDisplacementsSame_v2(displacements1sIncrement, uInitialFreeDOFs_state1.CopyToArray()));
            Assert.True(NRNLAnalyzerDevelopTest_v2.AreDisplacementsSame_v2(displacements2ndIncrement, uInitialFreeDOFs_state2.CopyToArray()));
            Assert.True(BondSlipTest.AreDisplacementsSame_v2(consCheck1, consCheck1Expected));
            Assert.True(NRNLAnalyzerDevelopTest_v2.AreDisplacementsSame_v2(stressesCheck3, stressesCheck3Expected));
            Assert.True(NRNLAnalyzerDevelopTest_v2.AreDisplacementsSame_v2(stressesCheck4, stressesCheck4Expected));
        }

        [Fact]
        public static void CheckOneRveParallel()
        {
            (int[] hexaPrint,int[] cohePrint, int[] shellPrint) = OneRveExample.Check_Graphene_rve_parallel();

            string results_file1 = "..\\..\\..\\RveTemplates\\Input\\RveGrShMultiple\\rve_no_1\\subdomainHexas.txt";
            string results_file2 = "..\\..\\..\\RveTemplates\\Input\\RveGrShMultiple\\rve_no_1\\subdomainCohesiveElements.txt";
            string results_file3 = "..\\..\\..\\RveTemplates\\Input\\RveGrShMultiple\\rve_no_1\\subdomainShellElements.txt";
            //string results_file4 = "..\\..\\..\\InputFiles\\MicroStructureAndIntegrationOneRveMultipleGrSh\\stressesCheck3.txt";
            //string results_file5 = "..\\..\\..\\InputFiles\\MicroStructureAndIntegrationOneRveMultipleGrSh\\stressesCheck4.txt";
            int[] hexaPrintExpected = PrintUtilities.ReadIntVector(results_file1);
            int[] cohePrintExpected = PrintUtilities.ReadIntVector(results_file2);
            int[] shellPrintExpected = PrintUtilities.ReadIntVector(results_file3);
            

            Assert.True(AreDisplacementsSame_v2(hexaPrint, hexaPrintExpected));
            Assert.True(AreDisplacementsSame_v2(cohePrint, cohePrintExpected));
            Assert.True(AreDisplacementsSame_v2(shellPrint, shellPrintExpected));
        }

        public static bool AreDisplacementsSame_v2(int[] expectedValues,
            int[] computedValues)
        {
            var comparer = new ValueComparer(1E-14);
            for (int i1 = 0; i1 < expectedValues.GetLength(0); i1++)
            {
                if (!comparer.AreEqual(expectedValues[i1], computedValues[i1]))
                {
                    return false;
                }
            }
            return true;
        }
    }
}
