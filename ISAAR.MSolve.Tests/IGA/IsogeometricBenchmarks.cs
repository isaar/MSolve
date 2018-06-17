using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using ISAAR.MSolve.Analyzers;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.IGA;
using ISAAR.MSolve.IGA.Entities;
using ISAAR.MSolve.IGA.Entities.Loads;
using ISAAR.MSolve.IGA.Readers;
using ISAAR.MSolve.Materials;
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
	    public void IsogeometricQuadraticCantilever2DWithDistributedLoad()
	    {
			// Model
			VectorExtensions.AssignTotalAffinityCount();
			Model model = new Model();
			ModelCreator modelCreator = new ModelCreator(model);
			string filename = "..\\..\\..\\IGA\\InputFiles\\Cantilever2D.txt";
			IsogeometricReader modelReader = new IsogeometricReader(modelCreator, filename);
			modelReader.CreateModelFromFile();

			// Forces and Boundary Conditions
		    Value verticalDistributedLoad = delegate (double x, double y, double z)
		    {
			    return new double[] {0,-100,0};
		    };
			model.PatchesDictionary[0].EdgesDictionary[1].LoadingConditions.Add(new NeumannBoundaryCondition(verticalDistributedLoad));

			// Boundary Conditions - Dirichlet
			foreach (ControlPoint controlPoint in model.PatchesDictionary[0].EdgesDictionary[0].ControlPointsDictionary.Values)
			{
				model.ControlPointsDictionary[controlPoint.ID].Constrains.Add(DOFType.X);
				model.ControlPointsDictionary[controlPoint.ID].Constrains.Add(DOFType.Y);
			}
			model.ConnectDataStructures();

			// Solvers
			var linearSystems = new Dictionary<int, ILinearSystem>();
			linearSystems[0] = new SkylineLinearSystem(0, model.PatchesDictionary[0].Forces);
			SolverSkyline solver = new SolverSkyline(linearSystems[0]);
			ProblemStructural provider = new ProblemStructural(model, linearSystems);
			LinearAnalyzer analyzer = new LinearAnalyzer(solver, linearSystems);
			StaticAnalyzer parentAnalyzer = new StaticAnalyzer(provider, analyzer, linearSystems);

			parentAnalyzer.BuildMatrices();
			parentAnalyzer.Initialize();
			parentAnalyzer.Solve();


			//Test for Load Vector
		    double[] loadVectorExpected =
		    {
			    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
			    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
			    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -83.3333333333334, 0, -166.666666666667, 0,
			    -250.000000000000, 0, -250.000000000000, 0, -166.666666666667, 0, -83.3333333333334
		    };

		    for (int i = 0; i < loadVectorExpected.Length; i++)
				Assert.Equal(loadVectorExpected[i], model.PatchesDictionary[0].Forces[i],6);

			//Test fro Displacement Vector
		    double[] displacementVectorExpected = new double[]
		    {
			    -0.0614562631781765, -0.0219170576325587, -0.0393637953297948, -0.00985891244520424, -0.0126687982396880,
			    -0.000943275549340603, 0.0126687982396881, -0.000943275549340565, 0.0393637953297948, -0.00985891244520429,
			    0.0614562631781765, -0.0219170576325587, -0.164630630046951, -0.143599796794095, -0.122284459518038,
			    -0.137572020341543, -0.0396482606150141, -0.129909984638882, 0.0396482606150142, -0.129909984638882,
			    0.122284459518038, -0.137572020341543, 0.164630630046951, -0.143599796794095, -0.255126649525807,
			    -0.367249906501344, -0.188034807331162, -0.360563789448209, -0.0622102833016037, -0.354689958832955,
			    0.0622102833016041, -0.354689958832955, 0.188034807331162, -0.360563789448210, 0.255126649525807,
			    -0.367249906501344, -0.329973104647066, -0.672298362080536, -0.245006560671741, -0.667431304467007,
			    -0.0808025488806688, -0.662174063778831, 0.0808025488806692, -0.662174063778831, 0.245006560671741,
			    -0.667431304467006, 0.329973104647066, -0.672298362080536, -0.389985650923780, -1.04507710416196,
			    -0.289753557243170, -1.04105642269944, -0.0959521766831420, -1.03721560192611, 0.0959521766831417,
			    -1.03721560192610, 0.289753557243169, -1.04105642269944, 0.389985650923779, -1.04507710416196,
			    -0.435055143301642, -1.47030303359559, -0.323516436156270, -1.46754420136883, -0.107046279714281,
			    -1.46472910091973, 0.107046279714280, -1.46472910091973, 0.323516436156269, -1.46754420136883, 0.435055143301642,
			    -1.47030303359559, -0.464895102729005, -1.93322992417900, -0.346196414935151, -1.93146777511281,
			    -0.114797364187326, -1.92958905082341, 0.114797364187326, -1.92958905082341, 0.346196414935151,
			    -1.93146777511282, 0.464895102729004, -1.93322992417900, -0.480254160241311, -2.41767489777727,
			    -0.357018220686884, -2.41725889060713, -0.117854053323470, -2.41770706308454, 0.117854053323469,
			    -2.41770706308454, 0.357018220686882, -2.41725889060713, 0.480254160241309, -2.41767489777726,
			    -0.481234435380421, -2.66710696029277, -0.357129997781511, -2.66655837235305, -0.118099002572320,
			    -2.66369484882849, 0.118099002572319, -2.66369484882849, 0.357129997781510, -2.66655837235305, 0.481234435380421,
			    -2.66710696029277
		    };
			for (int i = 0; i < displacementVectorExpected.Length; i++)
				Assert.Equal(displacementVectorExpected[i], linearSystems[0].Solution[i], 6);
		}

        //[Fact]
        //public void TSplinesShellsBenchmark()
        //{
        //    VectorExtensions.AssignTotalAffinityCount();
        //    Model model = new Model();
        //    string filename = "..\\..\\..\\IGA\\InputFiles\\tspline.iga";
        //    IGAFileReader modelReader = new IGAFileReader(model, filename);
        //    modelReader.CreateTSplineShellsModelFromFile();

        //    model.PatchesDictionary[0].Material = new ElasticMaterial2D
        //    {
        //        PoissonRatio = 0.3,
        //        YoungModulus = 10e6,
        //        StressState = StressStates.PlaneStress
        //    };
        //    model.PatchesDictionary[0].Thickness = 0.1;

        //    for (int i = 0; i < 100; i++)
        //    {
        //        model.ControlPointsDictionary[i].Constrains.Add(DOFType.X);
        //        model.ControlPointsDictionary[i].Constrains.Add(DOFType.Y);
        //        model.ControlPointsDictionary[i].Constrains.Add(DOFType.Z);
        //    }

        //    for (int i = 0; i < model.ControlPoints.Count-100; i++)
        //    {
        //        model.Loads.Add(new Load() { Amount = -1000, ControlPoint = model.ControlPointsDictionary[i], DOF = DOFType.Y });
        //    }


        //    //model.Loads.Add(new Load() { Amount = -10000, ControlPoint = model.ControlPointsDictionary[2], DOF = DOFType.Y });
        //    model.ConnectDataStructures();

        //    var linearSystems = new Dictionary<int, ILinearSystem>();
        //    linearSystems[0] = new SkylineLinearSystem(0, model.PatchesDictionary[0].Forces);
        //    SolverSkyline solver = new SolverSkyline(linearSystems[0]);
        //    ProblemStructural provider = new ProblemStructural(model, linearSystems);
        //    LinearAnalyzer analyzer = new LinearAnalyzer(solver, linearSystems);
        //    StaticAnalyzer parentAnalyzer = new StaticAnalyzer(provider, analyzer, linearSystems);

        //    parentAnalyzer.BuildMatrices();
        //    parentAnalyzer.Initialize();
        //    parentAnalyzer.Solve();
        //}

    }
}
