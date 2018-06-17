using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using ISAAR.MSolve.Analyzers;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.IGA;
using ISAAR.MSolve.IGA.Entities;
using ISAAR.MSolve.IGA.Entities.Loads;
using ISAAR.MSolve.IGA.Problems.Structural.Constitutive;
using ISAAR.MSolve.IGA.Readers;
using ISAAR.MSolve.Numerical.LinearAlgebra;
using ISAAR.MSolve.Problems;
using ISAAR.MSolve.Solvers.Interfaces;
using ISAAR.MSolve.Solvers.Skyline;
using Xunit;

namespace ISAAR.MSolve.Tests.IGA
{
    public class IsogeometricBenchmarks
    {
		//[Fact]
		//public void Qube64()
		//{
		//	// Model
		//	VectorExtensions.AssignTotalAffinityCount();
		//	Model model = new Model();
		//	ModelCreator modelCreator = new ModelCreator(model);
		//	string filename = "..\\..\\..\\IGA\\InputFiles\\Qube4096.txt";
		//	IsogeometricReader modelReader = new IsogeometricReader(modelCreator, filename);
		//	modelReader.CreateModelFromFile();

		//	// Forces and Boundary Conditions
		//	// Loading Conditions - Pressure

		//	model.PatchesDictionary[0].FacesDictionary[3].LoadingConditions.Add(new PressureBoundaryCondition(100));

		//	// Boundary Conditions - Dirichlet
		//	foreach (ControlPoint controlPoint in model.PatchesDictionary[0].FacesDictionary[2].ControlPointsDictionary.Values)
		//	{
		//		model.ControlPointsDictionary[controlPoint.ID].Constrains.Add(DOFType.X);
		//		model.ControlPointsDictionary[controlPoint.ID].Constrains.Add(DOFType.Y);
		//		model.ControlPointsDictionary[controlPoint.ID].Constrains.Add(DOFType.Z);
		//	}
		//	model.ConnectDataStructures();

		//	System.IO.File.WriteAllLines("..\\..\\..\\IGA\\InputFiles\\Qube4096Forces.txt",
		//		model.PatchesDictionary[0].Forces.Select(f=>f.ToString()));
			
		//}

		[Fact]
        public void IsogeometricQuadraticCantilever2D()
        {
            // Model
            VectorExtensions.AssignTotalAffinityCount();
            Model model = new Model();
            ModelCreator modelCreator = new ModelCreator(model);
            string filename = "..\\..\\..\\IGA\\InputFiles\\Cantilever2D.txt";
            IsogeometricReader modelReader = new IsogeometricReader(modelCreator, filename);
            modelReader.CreateModelFromFile();            

            // Forces and Boundary Conditions
            foreach (ControlPoint controlPoint in model.PatchesDictionary[0].EdgesDictionary[1].ControlPointsDictionary.Values)
                model.Loads.Add(new Load() { Amount = -100, ControlPoint = model.ControlPoints[controlPoint.ID], DOF = DOFType.Y });

            // Boundary Conditions - Dirichlet
            foreach (ControlPoint controlPoint in model.PatchesDictionary[0].EdgesDictionary[0].ControlPointsDictionary.Values)
            {
                model.ControlPointsDictionary[controlPoint.ID].Constrains.Add(DOFType.X);
                model.ControlPointsDictionary[controlPoint.ID].Constrains.Add(DOFType.Y);
            }
            model.ConnectDataStructures();

            // Solvers
            var linearSystems = new Dictionary<int, ILinearSystem>(); //I think this should be done automatically
            linearSystems[0] = new SkylineLinearSystem(0, model.PatchesDictionary[0].Forces);
            SolverSkyline solver = new SolverSkyline(linearSystems[0]);
            ProblemStructural provider = new ProblemStructural(model, linearSystems);
            LinearAnalyzer analyzer = new LinearAnalyzer(solver, linearSystems);
            StaticAnalyzer parentAnalyzer = new StaticAnalyzer(provider, analyzer, linearSystems);

            parentAnalyzer.BuildMatrices();
            parentAnalyzer.Initialize();
            parentAnalyzer.Solve();


            //Test for Load Vector
            int[] loadVectorExpected = new int[108];
            loadVectorExpected[110 - 13] = -100;
            loadVectorExpected[112 - 13] = -100;
            loadVectorExpected[114 - 13] = -100;
            loadVectorExpected[116 - 13] = -100;
            loadVectorExpected[118 - 13] = -100;
            loadVectorExpected[120 - 13] = -100;

            for (int i = 0; i < loadVectorExpected.Length; i++)
                Assert.Equal(loadVectorExpected[i], model.PatchesDictionary[0].Forces[i]);

            //Test fro Displacement Vector
            double[] displacementVectorExpected = new double[]
            {
                -0.0368738351302207, -0.0131499592237545, -0.0236187564608231,
                -0.00591458123276355, -0.00760062703427646, -0.000566403168697526, 0.00760062703427637,
                -0.000566403168697477, 0.0236187564608231, -0.00591458123276345, 0.0368738351302207,
                -0.0131499592237545, -0.0987784914203237, -0.0861605620323561, -0.0733694959613918, -0.0825449187386709,
                -0.0237898071790881, -0.0779453691704672, 0.0237898071790880, -0.0779453691704673, 0.0733694959613916,
                -0.0825449187386711, 0.0987784914203237, -0.0861605620323562, -0.153074952778768, -0.220348665003477,
                -0.112824242793523, -0.216334418527862, -0.0373252643977994, -0.212814901275071, 0.0373252643977993,
                -0.212814901275071, 0.112824242793523, -0.216334418527862, 0.153074952778767, -0.220348665003477,
                -0.197987950979100, -0.403380316332308, -0.146994646367387, -0.400466656395154, -0.0484809128603281,
                -0.397304037977320, 0.0484809128603277, -0.397304037977320, 0.146994646367387, -0.400466656395153,
                0.197987950979101, -0.403380316332307, -0.233979100719325, -0.627050378576386, -0.173875915391434,
                -0.624621548530984, -0.0575789187734497, -0.622325573362235, 0.0575789187734493, -0.622325573362235,
                0.173875915391433, -0.624621548530985, 0.233979100719325, -0.627050378576387, -0.261062696370552,
                -0.882149456584904, -0.194054136046524, -0.880534402771693, -0.0641987082249154, -0.878857804164279,
                0.0641987082249148, -0.878857804164278, 0.194054136046524, -0.880534402771692, 0.261062696370551,
                -0.882149456584903, -0.278878750048815, -1.16007242947546, -0.207837054153937, -1.15891585124500,
                -0.0689580379467459, -1.15768229264197, 0.0689580379467447, -1.15768229264197, 0.207837054153936,
                -1.15891585124500, 0.278878750048814, -1.16007242947546, -0.288294648281518, -1.45007762926002,
                -0.213983781845837, -1.45008375093483, -0.0704137097578954, -1.45093051606850, 0.0704137097578944,
                -1.45093051606850, 0.213983781845836, -1.45008375093483, 0.288294648281518, -1.45007762926002,
                -0.289593315112231, -1.60204696081608, -0.214138050751406, -1.60097881234833, -0.0706978642320891,
                -1.59760921075339, 0.0706978642320882, -1.59760921075338, 0.214138050751405, -1.60097881234833,
                0.289593315112230, -1.60204696081608
            };

            for (int i = 0; i < displacementVectorExpected.Length; i++)
                Assert.Equal(displacementVectorExpected[i], linearSystems[0].Solution[i], 6);

        }

        [Fact]
        public void TSplinesShellsBenchmark()
        {
            VectorExtensions.AssignTotalAffinityCount();
            Model model = new Model();
            string filename = "..\\..\\..\\IGA\\InputFiles\\tspline.iga";
            IGAFileReader modelReader = new IGAFileReader(model, filename);
            modelReader.CreateTSplineShellsModelFromFile();

            model.PatchesDictionary[0].Material = new ElasticMaterial2D
            {
                PoissonRatio = 0.3,
                YoungModulus = 10e6,
                StressState = StressStates.PlaneStress
            };
            model.PatchesDictionary[0].Thickness = 0.1;

            for (int i = 0; i < 100; i++)
            {
                model.ControlPointsDictionary[i].Constrains.Add(DOFType.X);
                model.ControlPointsDictionary[i].Constrains.Add(DOFType.Y);
                model.ControlPointsDictionary[i].Constrains.Add(DOFType.Z);
            }

            for (int i = 0; i < model.ControlPoints.Count-100; i++)
            {
                model.Loads.Add(new Load() { Amount = -1000, ControlPoint = model.ControlPointsDictionary[i], DOF = DOFType.Y });
            }


            //model.Loads.Add(new Load() { Amount = -10000, ControlPoint = model.ControlPointsDictionary[2], DOF = DOFType.Y });
            model.ConnectDataStructures();

            var linearSystems = new Dictionary<int, ILinearSystem>();
            linearSystems[0] = new SkylineLinearSystem(0, model.PatchesDictionary[0].Forces);
            SolverSkyline solver = new SolverSkyline(linearSystems[0]);
            ProblemStructural provider = new ProblemStructural(model, linearSystems);
            LinearAnalyzer analyzer = new LinearAnalyzer(solver, linearSystems);
            StaticAnalyzer parentAnalyzer = new StaticAnalyzer(provider, analyzer, linearSystems);

            parentAnalyzer.BuildMatrices();
            parentAnalyzer.Initialize();
            parentAnalyzer.Solve();
        }

    }
}
