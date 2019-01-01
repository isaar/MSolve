using System;
using System.Collections.Generic;
using ISAAR.MSolve.Analyzers;
using ISAAR.MSolve.Discretization;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.IGA.Elements;
using ISAAR.MSolve.IGA.Entities;
using ISAAR.MSolve.IGA.Entities.Loads;
using ISAAR.MSolve.IGA.Postprocessing;
using ISAAR.MSolve.IGA.Readers;
using ISAAR.MSolve.Materials;
using ISAAR.MSolve.Numerical.LinearAlgebra;
using ISAAR.MSolve.Problems;
using ISAAR.MSolve.Solvers.Direct;
using ISAAR.MSolve.Solvers.Interfaces;
using ISAAR.MSolve.Solvers.Ordering;
using ISAAR.MSolve.Solvers.Ordering.Reordering;
using ISAAR.MSolve.Solvers.Skyline;
using MathNet.Numerics.Data.Matlab;
using MathNet.Numerics.LinearAlgebra;
using Xunit;
using VectorExtensions = ISAAR.MSolve.Numerical.LinearAlgebra.VectorExtensions;

namespace ISAAR.MSolve.IGA.Tests
{
	public class IsogeometricBenchmarks
	{
		[Fact]
		public void IsogeometricQuadraticCantilever2D()
		{
			// Model
			VectorExtensions.AssignTotalAffinityCount();
			Model model = new Model();
			ModelCreator modelCreator = new ModelCreator(model);
			string filename = "..\\..\\..\\InputFiles\\Cantilever2D.txt";
			IsogeometricReader modelReader = new IsogeometricReader(modelCreator, filename);
			modelReader.CreateModelFromFile();

			// Forces and Boundary Conditions
			foreach (ControlPoint controlPoint in model.PatchesDictionary[0].EdgesDictionary[1].ControlPointsDictionary
				.Values)
				model.Loads.Add(new Load()
				{
					Amount = -100,
					ControlPoint = model.ControlPoints[controlPoint.ID],
					DOF = DOFType.Y
				});

			// Boundary Conditions - Dirichlet
			foreach (ControlPoint controlPoint in model.PatchesDictionary[0].EdgesDictionary[0].ControlPointsDictionary
				.Values)
			{
				model.ControlPointsDictionary[controlPoint.ID].Constrains.Add(new Constraint() {DOF = DOFType.X});
				model.ControlPointsDictionary[controlPoint.ID].Constrains.Add(new Constraint() { DOF = DOFType.Y });
			}

			var solverBuilder = new DenseMatrixSolver.Builder();
			solverBuilder.DofOrderer = new DofOrderer(
				new NodeMajorDofOrderingStrategy(), new NullReordering());
			ISolver_v2 solver = solverBuilder.BuildSolver(model);

			// Structural problem provider
			var provider = new ProblemStructural_v2(model, solver);

			// Linear static analysis
			var childAnalyzer = new LinearAnalyzer_v2(solver);
			var parentAnalyzer = new StaticAnalyzer_v2(model, solver, provider, childAnalyzer);

			// Run the analysis
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
				Assert.Equal(displacementVectorExpected[i], solver.LinearSystems[0].Solution[i], 6);
		}

		[Fact]
		public void IsogeometricQuadraticCantilever2DWithDistributedLoad()
		{
			// Model
			VectorExtensions.AssignTotalAffinityCount();
			Model model = new Model();
			ModelCreator modelCreator = new ModelCreator(model);
			var filename = "Cantilever2D";
			string filepath =$"..\\..\\..\\InputFiles\\{filename}.txt";
			IsogeometricReader modelReader = new IsogeometricReader(modelCreator, filepath);
			modelReader.CreateModelFromFile();

			// Forces and Boundary Conditions
			Value verticalDistributedLoad = delegate(double x, double y, double z)
			{
				return new double[] {0, -100, 0};
			};
			model.PatchesDictionary[0].EdgesDictionary[1].LoadingConditions
				.Add(new NeumannBoundaryCondition(verticalDistributedLoad));

			// Boundary Conditions - Dirichlet
			foreach (ControlPoint controlPoint in model.PatchesDictionary[0].EdgesDictionary[0].ControlPointsDictionary
				.Values)
			{
				model.ControlPointsDictionary[controlPoint.ID].Constrains.Add(new Constraint() {DOF = DOFType.X});
				model.ControlPointsDictionary[controlPoint.ID].Constrains.Add(new Constraint() { DOF = DOFType.Y });
			}

			var solverBuilder = new SuiteSparseSolver.Builder();
			solverBuilder.DofOrderer = new DofOrderer(
				new NodeMajorDofOrderingStrategy(), new NullReordering());
			ISolver_v2 solver = solverBuilder.BuildSolver(model);

			// Structural problem provider
			var provider = new ProblemStructural_v2(model, solver);

			// Linear static analysis
			var childAnalyzer = new LinearAnalyzer_v2(solver);
			var parentAnalyzer = new StaticAnalyzer_v2(model, solver, provider, childAnalyzer);

			// Run the analysis
			parentAnalyzer.Initialize();
			parentAnalyzer.Solve();

			var paraviewOutput = new ParaviewNurbs2D(model,solver.LinearSystems[0].Solution, filename);
			paraviewOutput.CreateParaview2DFile();


			//Test for Load Vector
			double[] loadVectorExpected =
			{
				0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
				0, 0, 0,
				0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
				0, 0, 0,
				0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -83.3333333333334, 0, -166.666666666667,
				0,
				-250.000000000000, 0, -250.000000000000, 0, -166.666666666667, 0, -83.3333333333334
			};

			for (int i = 0; i < loadVectorExpected.Length; i++)
				Assert.Equal(loadVectorExpected[i], model.PatchesDictionary[0].Forces[i], 6);

			//Test fro Displacement Vector
			double[] displacementVectorExpected = new double[]
			{
				-0.0614562631781765, -0.0219170576325587, -0.0393637953297948, -0.00985891244520424,
				-0.0126687982396880,
				-0.000943275549340603, 0.0126687982396881, -0.000943275549340565, 0.0393637953297948,
				-0.00985891244520429,
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
				-1.46472910091973, 0.107046279714280, -1.46472910091973, 0.323516436156269, -1.46754420136883,
				0.435055143301642,
				-1.47030303359559, -0.464895102729005, -1.93322992417900, -0.346196414935151, -1.93146777511281,
				-0.114797364187326, -1.92958905082341, 0.114797364187326, -1.92958905082341, 0.346196414935151,
				-1.93146777511282, 0.464895102729004, -1.93322992417900, -0.480254160241311, -2.41767489777727,
				-0.357018220686884, -2.41725889060713, -0.117854053323470, -2.41770706308454, 0.117854053323469,
				-2.41770706308454, 0.357018220686882, -2.41725889060713, 0.480254160241309, -2.41767489777726,
				-0.481234435380421, -2.66710696029277, -0.357129997781511, -2.66655837235305, -0.118099002572320,
				-2.66369484882849, 0.118099002572319, -2.66369484882849, 0.357129997781510, -2.66655837235305,
				0.481234435380421,
				-2.66710696029277
			};
			for (int i = 0; i < displacementVectorExpected.Length; i++)
				Assert.Equal(displacementVectorExpected[i], solver.LinearSystems[0].Solution[i], 6);
		}

		[Fact]
		public void IsogeometricPlateTension()
		{
			// Model
			VectorExtensions.AssignTotalAffinityCount();
			Model model = new Model();
			ModelCreator modelCreator = new ModelCreator(model);
			string filename = "..\\..\\..\\InputFiles\\PlateTension.txt";
			IsogeometricReader modelReader = new IsogeometricReader(modelCreator, filename);
			modelReader.CreateModelFromFile();

			// Forces and Boundary Conditions
			Value horizontalDistributedLoad = delegate(double x, double y, double z)
			{
				return new double[] {100, 0, 0};
			};
			model.PatchesDictionary[0].EdgesDictionary[1].LoadingConditions
				.Add(new NeumannBoundaryCondition(horizontalDistributedLoad));

			// Boundary Conditions - Dirichlet
			foreach (ControlPoint controlPoint in model.PatchesDictionary[0].EdgesDictionary[0].ControlPointsDictionary
				.Values)
			{
				model.ControlPointsDictionary[controlPoint.ID].Constrains.Add(new Constraint() {DOF = DOFType.X});
			}

			model.ControlPointsDictionary[0].Constrains.Add(new Constraint() { DOF = DOFType.Y });

			var solverBuilder = new DenseMatrixSolver.Builder();
			solverBuilder.DofOrderer = new DofOrderer(
				new NodeMajorDofOrderingStrategy(), new NullReordering());
			ISolver_v2 solver = solverBuilder.BuildSolver(model);

			// Structural problem provider
			var provider = new ProblemStructural_v2(model, solver);

			// Linear static analysis
			var childAnalyzer = new LinearAnalyzer_v2(solver);
			var parentAnalyzer = new StaticAnalyzer_v2(model, solver, provider, childAnalyzer);

			// Run the analysis
			parentAnalyzer.Initialize();
			parentAnalyzer.Solve();

			double[] forceVectorExpected =
			{
				0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
				0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
				0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
				0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
				0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
				0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
				2.813476563, 0, 5.638671875, 0, 8.510742188, 0, 11.50390625, 0, 11.796875, 0, 12.20703125, 0, 12.734375,
				0, 13.37890625, 0, 10.47949219, 0, 7.232421875, 0, 3.704101563, 0,
			};
			for (int i = 0; i < forceVectorExpected.Length; i++)
				Assert.Equal(forceVectorExpected[i], model.PatchesDictionary[0].Forces[i], 8);

			#region expectedDisplacement

			double[] displacementVectorExpected =
			{
				-0.000146256, -0.000438751, -0.000882074, -0.001334534, -0.001800705, -0.002285162, -0.00279246,
				-0.003327199, -0.003704999, -0.003899999, 0.000341248, -7.25E-09, 0.000341253, -0.00014625, 0.000341246,
				-0.000438757, 0.000341253, -0.00088207, 0.000341247, -0.001334537, 0.000341252, -0.001800705,
				0.000341249, -0.002285157, 0.000341247, -0.002792472, 0.000341263, -0.003327173, 0.000341224,
				-0.003705018, 0.000341261, -0.003900014, 0.001023752, -3.43E-13, 0.001023744, -0.000146257, 0.001023757,
				-0.000438746, 0.001023742, -0.000882079, 0.001023758, -0.001334528, 0.001023743, -0.001800711,
				0.001023754, -0.002285159, 0.001023755, -0.002792456, 0.001023724, -0.003327216, 0.001023809,
				-0.00370498, 0.001023732, -0.003899995, 0.00205816, -6.57E-09, 0.002058168, -0.00014625, 0.002058153,
				-0.000438762, 0.002058174, -0.000882065, 0.002058151, -0.001334544, 0.002058177, -0.001800696,
				0.002058155, -0.002285166, 0.002058162, -0.002792467, 0.002058201, -0.003327172, 0.002058075,
				-0.00370503, 0.002058194, -0.003900003, 0.003113908, 5.09E-10, 0.003113896, -0.000146258, 0.003113922,
				-0.000438741, 0.003113886, -0.000882088, 0.003113929, -0.001334515, 0.00311388, -0.001800728,
				0.00311393, -0.002285137, 0.003113904, -0.002792481, 0.003113852, -0.003327193, 0.00311407,
				-0.003704984, 0.003113858, -0.003900043, 0.004201632, -1.28E-08, 0.004201651, -0.00014625, 0.004201609,
				-0.000438773, 0.004201673, -0.000882051, 0.004201595, -0.00133457, 0.004201694, -0.001800657,
				0.00420159, -0.00228522, 0.004201668, -0.002792395, 0.004201722, -0.003327265, 0.004201363,
				-0.003704993, 0.004201716, -0.003899863, 0.00533204, 9.73E-09, 0.005332006, -0.000146258, 0.005332079,
				-0.000438732, 0.005331968, -0.000882113, 0.005332111, -0.001334471, 0.005331923, -0.001800805,
				0.005332149, -0.00228501, 0.005331957, -0.002792662, 0.005331965, -0.003326937, 0.005332439,
				-0.003705134, 0.005331965, -0.003900417, 0.006515714, -3.75E-08, 0.006515781, -0.00014626, 0.006515657,
				-0.000438787, 0.006515838, -0.000882016, 0.006515598, -0.001334655, 0.006515925, -0.001800499,
				0.006515523, -0.002285491, 0.006515942, -0.002791946, 0.006515783, -0.003327922, 0.006515237,
				-0.003704656, 0.006515668, -0.003898999, 0.007763498, 5.72E-08, 0.007763343, -0.000146237, 0.007763593,
				-0.000438734, 0.00776324, -0.00088218, 0.007763706, -0.001334301, 0.007763035, -0.001801172,
				0.007763936, -0.002284351, 0.007762843, -0.002793734, 0.007764044, -0.003325117, 0.007763485,
				-0.003706045, 0.007764022, -0.003902662, 0.008644866, -1.29E-07, 0.00864517, -0.000146322, 0.008644803,
				-0.000438748, 0.008645196, -0.000881942, 0.00864469, -0.001334917, 0.008645477, -0.001799884,
				0.008644261, -0.00228676, 0.008646, -0.00278979, 0.008643379, -0.003331508, 0.008647443, -0.003702447,
				0.008641489, -0.003894147, 0.00910031, 1.35E-07, 0.009099715, -0.000146173, 0.009100131, -0.000438822,
				0.009099933, -0.000882057, 0.009100099, -0.001334479, 0.009099651, -0.001800931, 0.009101019,
				-0.002284642, 0.009097205, -0.00279357, 0.009107416, -0.003325298, 0.009086905, -0.003706366,
				0.009112171, -0.003906529
			};
			for (int i = 0; i < displacementVectorExpected.Length; i++)
				Assert.Equal(displacementVectorExpected[i], solver.LinearSystems[0].Solution[i], 4);

			#endregion
		}
		
		[Fact]
		public void IsogeometricPlaneStrainMixedBC()
		{
			// Model
			VectorExtensions.AssignTotalAffinityCount();
			Model model = new Model();
			ModelCreator modelCreator = new ModelCreator(model);
			string filename = "..\\..\\..\\InputFiles\\SquareMixedBC.txt";
			IsogeometricReader modelReader = new IsogeometricReader(modelCreator, filename);
			modelReader.CreateModelFromFile();

			// Forces and Boundary Conditions
			//model.PatchesDictionary[0].EdgesDictionary[1].LoadingConditions.Add(new PressureBoundaryCondition(0.3));
			//model.PatchesDictionary[0].EdgesDictionary[3].LoadingConditions.Add(new PressureBoundaryCondition(0.3));
			Value horizontalDistributedLoad = delegate(double x, double y, double z)
			{
				return new double[] {-0.30, 0, 0};
			};
			model.PatchesDictionary[0].EdgesDictionary[1].LoadingConditions
				.Add(new NeumannBoundaryCondition(horizontalDistributedLoad));

			Value verticalDistributedLoad = delegate(double x, double y, double z)
			{
				return new double[] {0, -0.3, 0};
			};
			model.PatchesDictionary[0].EdgesDictionary[3].LoadingConditions
				.Add(new NeumannBoundaryCondition(verticalDistributedLoad));

			// Boundary Conditions - Dirichlet
			foreach (ControlPoint controlPoint in model.PatchesDictionary[0].EdgesDictionary[0].ControlPointsDictionary
				.Values)
			{
				model.ControlPointsDictionary[controlPoint.ID].Constrains.Add(new Constraint() {DOF = DOFType.X});
				model.ControlPointsDictionary[controlPoint.ID].Constrains.Add(new Constraint() { DOF = DOFType.Y });
			}

			foreach (ControlPoint controlPoint in model.PatchesDictionary[0].EdgesDictionary[2].ControlPointsDictionary
				.Values)
			{
				model.ControlPointsDictionary[controlPoint.ID].Constrains.Add(new Constraint() { DOF = DOFType.X });
				model.ControlPointsDictionary[controlPoint.ID].Constrains.Add(new Constraint() { DOF = DOFType.Y });
			}
		

			// Solvers
			var solverBuilder = new SuiteSparseSolver.Builder();
			solverBuilder.DofOrderer = new DofOrderer(
				new NodeMajorDofOrderingStrategy(), new NullReordering());
			ISolver_v2 solver = solverBuilder.BuildSolver(model);

			// Structural problem provider
			var provider = new ProblemStructural_v2(model, solver);

			// Linear static analysis
			var childAnalyzer = new LinearAnalyzer_v2(solver);
			var parentAnalyzer = new StaticAnalyzer_v2(model, solver, provider, childAnalyzer);

			// Run the analysis
			parentAnalyzer.Initialize();
			parentAnalyzer.Solve();

			var forceVectorExpected = new double[]
			{
				0,
				0,
				0,
				0,
				0,
				0,
				0,
				0,
				0,
				0,
				0,
				0,
				0,
				0,
				0,
				0,
				0,
				0,
				0,
				0,
				0,
				-0.0166666666666667,
				0,
				0,
				0,
				0,
				0,
				0,
				0,
				0,
				0,
				0,
				0,
				0,
				0,
				0,
				0,
				0,
				0,
				0,
				0,
				0,
				0,
				-0.0250000000000000,
				0,
				0,
				0,
				0,
				0,
				0,
				0,
				0,
				0,
				0,
				0,
				0,
				0,
				0,
				0,
				0,
				0,
				0,
				0,
				0,
				0,
				-0.0333333333333333,
				0,
				0,
				0,
				0,
				0,
				0,
				0,
				0,
				0,
				0,
				0,
				0,
				0,
				0,
				0,
				0,
				0,
				0,
				0,
				0,
				0,
				-0.0333333333333333,
				0,
				0,
				0,
				0,
				0,
				0,
				0,
				0,
				0,
				0,
				0,
				0,
				0,
				0,
				0,
				0,
				0,
				0,
				0,
				0,
				0,
				-0.0333333333333333,
				0,
				0,
				0,
				0,
				0,
				0,
				0,
				0,
				0,
				0,
				0,
				0,
				0,
				0,
				0,
				0,
				0,
				0,
				0,
				0,
				0,
				-0.0333333333333333,
				0,
				0,
				0,
				0,
				0,
				0,
				0,
				0,
				0,
				0,
				0,
				0,
				0,
				0,
				0,
				0,
				0,
				0,
				0,
				0,
				0,
				-0.0333333333333333,
				0,
				0,
				0,
				0,
				0,
				0,
				0,
				0,
				0,
				0,
				0,
				0,
				0,
				0,
				0,
				0,
				0,
				0,
				0,
				0,
				0,
				-0.0333333333333333,
				0,
				0,
				0,
				0,
				0,
				0,
				0,
				0,
				0,
				0,
				0,
				0,
				0,
				0,
				0,
				0,
				0,
				0,
				0,
				0,
				0,
				-0.0250000000000000,
				0,
				0,
				0,
				0,
				0,
				0,
				0,
				0,
				0,
				0,
				0,
				0,
				0,
				0,
				0,
				0,
				0,
				0,
				0,
				0,
				0,
				-0.0166666666666667,
				-0.0166666666666667,
				0,
				-0.0250000000000000,
				0,
				-0.0333333333333333,
				0,
				-0.0333333333333333,
				0,
				-0.0333333333333333,
				0,
				-0.0333333333333333,
				0,
				-0.0333333333333333,
				0,
				-0.0333333333333333,
				0,
				-0.0250000000000000,
				0,
				-0.0166666666666667,
				0,
				-0.00833333333333334,
				-0.00833333333333334,
			};

			for (int i = 0; i < forceVectorExpected.Length; i++)
				Assert.Equal(forceVectorExpected[i], model.PatchesDictionary[0].Forces[i], 7);

			#region expectedDisplacement

			var displacementVectorExpected = new double[]
			{
				-0.00314820363909275,
				-0.00314820363909275,
				-0.00379395242120674,
				-0.00474046877963962,
				-0.00527896515562401,
				-0.00627080476725420,
				-0.00611202232103848,
				-0.00732618402014196,
				-0.00671664390713507,
				-0.00830517563434155,
				-0.00750778852140214,
				-0.00921569343748655,
				-0.00711656646035208,
				-0.0104552048087267,
				-0.00857783159759678,
				-0.0115132239445180,
				-0.00447714927679477,
				-0.0142326173251668,
				-0.00896507378072984,
				-0.0172326398522936,
				0.0165696344045362,
				-0.0344920417612453,
				-0.00474046877963962,
				-0.00379395242120674,
				-0.0106163667467970,
				-0.0106163667467970,
				-0.0145901385325659,
				-0.0159579381873314,
				-0.0175798662746167,
				-0.0200752820397931,
				-0.0197642416239734,
				-0.0232933594058257,
				-0.0209594659687415,
				-0.0269739141344681,
				-0.0221701240269120,
				-0.0304381009757610,
				-0.0204015380221106,
				-0.0372162332314755,
				-0.0190262809325074,
				-0.0451261112082839,
				0.000524263735174436,
				-0.0646122725733136,
				0.0134058371171140,
				-0.0674237539857168,
				-0.00627080476725419,
				-0.00527896515562400,
				-0.0159579381873313,
				-0.0145901385325660,
				-0.0255583394854834,
				-0.0255583394854834,
				-0.0320458295432382,
				-0.0340345779839745,
				-0.0364394696445877,
				-0.0416128738023956,
				-0.0392520264843806,
				-0.0489264479827333,
				-0.0394128785675202,
				-0.0579962012559210,
				-0.0375262954335112,
				-0.0685615165413566,
				-0.0269408733978005,
				-0.0856583155480079,
				-0.0140932604372300,
				-0.0958837290337786,
				-0.00425322389965576,
				-0.102864805751424,
				-0.00732618402014194,
				-0.00611202232103848,
				-0.0200752820397931,
				-0.0175798662746167,
				-0.0340345779839744,
				-0.0320458295432382,
				-0.0442399103415714,
				-0.0442399103415714,
				-0.0511182027891641,
				-0.0551752972826913,
				-0.0547664019599041,
				-0.0662589972721233,
				-0.0554544754480729,
				-0.0781539846420991,
				-0.0514092081933032,
				-0.0928498000772380,
				-0.0427832850822321,
				-0.108645057746992,
				-0.0311333714176029,
				-0.120963904936585,
				-0.0251089003632628,
				-0.125275523366600,
				-0.00830517563434156,
				-0.00671664390713508,
				-0.0232933594058256,
				-0.0197642416239734,
				-0.0416128738023955,
				-0.0364394696445877,
				-0.0551752972826911,
				-0.0511182027891642,
				-0.0644659012290671,
				-0.0644659012290672,
				-0.0694887264699124,
				-0.0777115338314980,
				-0.0702085098844087,
				-0.0920154063898269,
				-0.0670415920146972,
				-0.107342294443253,
				-0.0590775341097851,
				-0.124438661501351,
				-0.0522405165453920,
				-0.134720354624928,
				-0.0486371118189779,
				-0.139962781550140,
				-0.00921569343748652,
				-0.00750778852140214,
				-0.0269739141344680,
				-0.0209594659687415,
				-0.0489264479827331,
				-0.0392520264843806,
				-0.0662589972721231,
				-0.0547664019599042,
				-0.0777115338314978,
				-0.0694887264699125,
				-0.0840291771115493,
				-0.0840291771115496,
				-0.0855834304606439,
				-0.0992380743950709,
				-0.0829894406989148,
				-0.115516212797889,
				-0.0776932540445211,
				-0.131918682864287,
				-0.0732614094641817,
				-0.142841837190617,
				-0.0714730797892666,
				-0.147708704990261,
				-0.0104552048087266,
				-0.00711656646035210,
				-0.0304381009757609,
				-0.0221701240269120,
				-0.0579962012559208,
				-0.0394128785675203,
				-0.0781539846420988,
				-0.0554544754480730,
				-0.0920154063898266,
				-0.0702085098844090,
				-0.0992380743950705,
				-0.0855834304606441,
				-0.101501531630957,
				-0.101501531630958,
				-0.0999370644850315,
				-0.118017394845313,
				-0.0962343442201912,
				-0.134878482285728,
				-0.0940899630707035,
				-0.145633866412365,
				-0.0936977503670338,
				-0.150902312951907,
				-0.0115132239445179,
				-0.00857783159759677,
				-0.0372162332314754,
				-0.0204015380221106,
				-0.0685615165413564,
				-0.0375262954335112,
				-0.0928498000772377,
				-0.0514092081933034,
				-0.107342294443253,
				-0.0670415920146973,
				-0.115516212797889,
				-0.0829894406989151,
				-0.118017394845313,
				-0.0999370644850319,
				-0.117164741332239,
				-0.117164741332239,
				-0.114732961045519,
				-0.134270598345563,
				-0.113764345865824,
				-0.145437001617946,
				-0.114021716983238,
				-0.150889885103904,
				-0.0142326173251668,
				-0.00447714927679481,
				-0.0451261112082838,
				-0.0190262809325074,
				-0.0856583155480077,
				-0.0269408733978006,
				-0.108645057746992,
				-0.0427832850822322,
				-0.124438661501351,
				-0.0590775341097853,
				-0.131918682864286,
				-0.0776932540445212,
				-0.134878482285727,
				-0.0962343442201915,
				-0.134270598345562,
				-0.114732961045519,
				-0.132671391727535,
				-0.132671391727535,
				-0.132190010013491,
				-0.144083210808925,
				-0.132499710574820,
				-0.149832975007404,
				-0.0172326398522934,
				-0.00896507378072984,
				-0.0646122725733136,
				0.000524263735174192,
				-0.0958837290337782,
				-0.0140932604372301,
				-0.120963904936585,
				-0.0311333714176032,
				-0.134720354624928,
				-0.0522405165453921,
				-0.142841837190616,
				-0.0732614094641819,
				-0.145633866412364,
				-0.0940899630707035,
				-0.145437001617945,
				-0.113764345865824,
				-0.144083210808926,
				-0.132190010013491,
				-0.143845673350073,
				-0.143845673350073,
				-0.143883483258064,
				-0.149631002634067,
				-0.0344920417612453,
				0.0165696344045361,
				-0.0674237539857166,
				0.0134058371171139,
				-0.102864805751424,
				-0.00425322389965601,
				-0.125275523366599,
				-0.0251089003632630,
				-0.139962781550140,
				-0.0486371118189782,
				-0.147708704990260,
				-0.0714730797892668,
				-0.150902312951907,
				-0.0936977503670341,
				-0.150889885103902,
				-0.114021716983238,
				-0.149832975007405,
				-0.132499710574819,
				-0.149631002634066,
				-0.143883483258064,
				-0.149670908292985,
				-0.149670908292985,
			};

			for (int i = 0; i < displacementVectorExpected.Length; i++)
				Assert.Equal(displacementVectorExpected[i], solver.LinearSystems[0].Solution[i], 7);

			#endregion
		}

		#region CantileverShellData
		private List<ControlPoint> ElementControlPoints()
		{
			return new List<ControlPoint>
			{
				new ControlPoint {ID = 0, X = 0.0, Y = 0.0, Z = 0.0, WeightFactor = 1.0},
				new ControlPoint {ID = 1, X = 0.0, Y = 0.5, Z = 0.0, WeightFactor = 1.0},
				new ControlPoint {ID = 2, X = 0.0, Y = 1.0, Z = 0.0, WeightFactor = 1.0},
				new ControlPoint {ID = 3, X = 16.66666667, Y = 0.0, Z = 0.0, WeightFactor = 1.0},
				new ControlPoint {ID = 4, X = 16.66666667, Y = 0.5, Z = 0.0, WeightFactor = 1.0},
				new ControlPoint {ID = 5, X = 16.66666667, Y = 1.0, Z = 0.0, WeightFactor = 1.0},
				new ControlPoint {ID = 6, X = 33.33333333, Y = 0.0, Z = 0.0, WeightFactor = 1.0},
				new ControlPoint {ID = 7, X = 33.33333333, Y = 0.5, Z = 0.0, WeightFactor = 1.0},
				new ControlPoint {ID = 8, X = 33.33333333, Y = 1.0, Z = 0.0, WeightFactor = 1.0},
				new ControlPoint {ID = 9, X = 50.0, Y = 0.0, Z = 0.0, WeightFactor = 1.0},
				new ControlPoint {ID = 10, X = 50.0, Y = 0.5, Z = 0.0, WeightFactor = 1.0},
				new ControlPoint {ID = 11, X = 50.0, Y = 1.0, Z = 0.0, WeightFactor = 1.0},
			};
		}

		private List<Knot> ElementKnot()
		{
			return new List<Knot>()
			{
				new Knot(){ID=0,Ksi=0.0,Heta=0.0,Zeta =0.0 },
				new Knot(){ID=1,Ksi=0.0,Heta=1.0,Zeta =0.0 },
				new Knot(){ID=2,Ksi=1.0,Heta=0.0,Zeta =0.0 },
				new Knot(){ID=3,Ksi=1.0,Heta=1.0,Zeta =0.0 }
			};
		}

		private Vector KnotValueVectorKsi()
		{
			return new Vector(new double[8]
			{
				0, 0, 0, 0, 1, 1, 1, 1
			});
		}

		private Vector KnotValueVectorHeta()
		{
			return new Vector(new double[6]
			{
				0, 0, 0, 1, 1, 1
			});
		}

		private NURBSKirchhoffLoveShellElement Element
		{
			get
			{
				var element = new NURBSKirchhoffLoveShellElement();
				var patch = new Patch();
				patch.Material = new ElasticMaterial2D(StressState2D.PlaneStrain)
				{
					YoungModulus = 100,
					PoissonRatio = 0.0
				};
				foreach (var controlPoint in ElementControlPoints())
					element.ControlPointsDictionary.Add(controlPoint.ID, controlPoint);
				foreach (var knot in ElementKnot())
					element.KnotsDictionary.Add(knot.ID, knot);
				element.ElementType = element;
				patch.Thickness = 1;
				patch.DegreeKsi = 3;
				patch.DegreeHeta = 2;
				patch.NumberOfControlPointsHeta = 3;
				patch.KnotValueVectorKsi = KnotValueVectorKsi();
				patch.KnotValueVectorHeta = KnotValueVectorHeta();
				element.Patch = patch;
				return element;
			}
		}
		#endregion
		
		//[Fact]
		public void IsogeometricHorseshoe3D()
		{
			VectorExtensions.AssignTotalAffinityCount();
			Model model = new Model();
			ModelCreator modelCreator = new ModelCreator(model);
			string filename = "..\\..\\..\\InputFiles\\Horseshoe.txt";
			IsogeometricReader modelReader = new IsogeometricReader(modelCreator, filename);
			modelReader.CreateModelFromFile();

			// Forces and Boundary Conditions
			Value verticalDistributedLoad = delegate (double x, double y, double z)
			{
				return new double[] { 0, 0, 0.1 };
			};
			model.PatchesDictionary[0].FacesDictionary[0].LoadingConditions.Add(new NeumannBoundaryCondition(verticalDistributedLoad));
			model.PatchesDictionary[0].FacesDictionary[1].LoadingConditions.Add(new NeumannBoundaryCondition(verticalDistributedLoad));
			model.PatchesDictionary[0].FacesDictionary[4].LoadingConditions.Add(new NeumannBoundaryCondition(verticalDistributedLoad));
			model.PatchesDictionary[0].FacesDictionary[5].LoadingConditions.Add(new NeumannBoundaryCondition(verticalDistributedLoad));


			// Boundary Conditions - Dirichlet
			foreach (ControlPoint controlPoint in model.PatchesDictionary[0].FacesDictionary[2].ControlPointsDictionary
				.Values)
			{
				model.ControlPointsDictionary[controlPoint.ID].Constrains.Add(new Constraint() {DOF = DOFType.X});
				model.ControlPointsDictionary[controlPoint.ID].Constrains.Add(new Constraint() { DOF = DOFType.Y });
				model.ControlPointsDictionary[controlPoint.ID].Constrains.Add(new Constraint() { DOF = DOFType.Z });
			}

			foreach (ControlPoint controlPoint in model.PatchesDictionary[0].FacesDictionary[3].ControlPointsDictionary
				.Values)
			{
				model.ControlPointsDictionary[controlPoint.ID].Constrains.Add(new Constraint() { DOF = DOFType.X });
				model.ControlPointsDictionary[controlPoint.ID].Constrains.Add(new Constraint() { DOF = DOFType.Y });
				model.ControlPointsDictionary[controlPoint.ID].Constrains.Add(new Constraint() { DOF = DOFType.Z });
			}

			// Solvers
			var solverBuilder = new SuiteSparseSolver.Builder();
			solverBuilder.DofOrderer = new DofOrderer(
				new NodeMajorDofOrderingStrategy(), new NullReordering());
			ISolver_v2 solver = solverBuilder.BuildSolver(model);

			// Structural problem provider
			var provider = new ProblemStructural_v2(model, solver);

			// Linear static analysis
			var childAnalyzer = new LinearAnalyzer_v2(solver);
			var parentAnalyzer = new StaticAnalyzer_v2(model, solver, provider, childAnalyzer);

			// Run the analysis
			parentAnalyzer.Initialize();
			parentAnalyzer.Solve();

			Matrix<double> forceVectorExpected = MatlabReader.Read<double>("..\\..\\..\\InputFiles\\Horseshoe.mat", "forceVector");

			for (int i = 0; i < forceVectorExpected.RowCount; i++)
				Assert.Equal(forceVectorExpected.At(i, 0), model.PatchesDictionary[0].Forces[i], 7);

			Matrix<double> displacementVectorExpected = MatlabReader.Read<double>("..\\..\\..\\InputFiles\\Horseshoe.mat", "displacementVector");

			for (int i = 0; i < displacementVectorExpected.RowCount; i++)
				Assert.Equal(displacementVectorExpected.At(i, 0), solver.LinearSystems[0].Solution[i], 7);

		}
		
		//[Fact]
		public void IsogeometricPlaneStrainRing()
		{
			VectorExtensions.AssignTotalAffinityCount();
			Model model = new Model();
			ModelCreator modelCreator = new ModelCreator(model);
			string filename = "..\\..\\..\\InputFiles\\SquareMixedBC.txt";
			IsogeometricReader modelReader = new IsogeometricReader(modelCreator, filename);
			modelReader.CreateModelFromFile();

			// Forces and Boundary Conditions
			//model.PatchesDictionary[0].EdgesDictionary[1].LoadingConditions.Add(new PressureBoundaryCondition(0.3));
			//model.PatchesDictionary[0].EdgesDictionary[3].LoadingConditions.Add(new PressureBoundaryCondition(0.3));
			Value horizontalDistributedLoad = delegate (double x, double y, double z)
			{
				return new double[] { -0.30, 0, 0 };
			};
			model.PatchesDictionary[0].EdgesDictionary[1].LoadingConditions
				.Add(new NeumannBoundaryCondition(horizontalDistributedLoad));

			Value verticalDistributedLoad = delegate (double x, double y, double z)
			{
				return new double[] { 0, -0.3, 0 };
			};
			model.PatchesDictionary[0].EdgesDictionary[3].LoadingConditions
				.Add(new NeumannBoundaryCondition(verticalDistributedLoad));

			// Boundary Conditions - Dirichlet
			foreach (ControlPoint controlPoint in model.PatchesDictionary[0].EdgesDictionary[0].ControlPointsDictionary
				.Values)
			{
				model.ControlPointsDictionary[controlPoint.ID].Constrains.Add(new Constraint() {DOF = DOFType.X});
				model.ControlPointsDictionary[controlPoint.ID].Constrains.Add(new Constraint() { DOF = DOFType.Y });
			}

			foreach (ControlPoint controlPoint in model.PatchesDictionary[0].EdgesDictionary[2].ControlPointsDictionary
				.Values)
			{
				model.ControlPointsDictionary[controlPoint.ID].Constrains.Add(new Constraint() { DOF = DOFType.X });
				model.ControlPointsDictionary[controlPoint.ID].Constrains.Add(new Constraint() { DOF = DOFType.Y });
			}

			// Solvers
			var solverBuilder = new SuiteSparseSolver.Builder();
			solverBuilder.DofOrderer = new DofOrderer(
				new NodeMajorDofOrderingStrategy(), new NullReordering());
			ISolver_v2 solver = solverBuilder.BuildSolver(model);

			// Structural problem provider
			var provider = new ProblemStructural_v2(model, solver);

			// Linear static analysis
			var childAnalyzer = new LinearAnalyzer_v2(solver);
			var parentAnalyzer = new StaticAnalyzer_v2(model, solver, provider, childAnalyzer);

			// Run the analysis
			parentAnalyzer.Initialize();
			parentAnalyzer.Solve();

			var forceVectorExpected = new double[]
			{
				0,
				0,
				0,
				0,
				0,
				0,
				0,
				0,
				0,
				0,
				0,
				0,
				0,
				0,
				0,
				0,
				0,
				0,
				0,
				0,
				0,
				-0.0166666666666667,
				0,
				0,
				0,
				0,
				0,
				0,
				0,
				0,
				0,
				0,
				0,
				0,
				0,
				0,
				0,
				0,
				0,
				0,
				0,
				0,
				0,
				-0.0250000000000000,
				0,
				0,
				0,
				0,
				0,
				0,
				0,
				0,
				0,
				0,
				0,
				0,
				0,
				0,
				0,
				0,
				0,
				0,
				0,
				0,
				0,
				-0.0333333333333333,
				0,
				0,
				0,
				0,
				0,
				0,
				0,
				0,
				0,
				0,
				0,
				0,
				0,
				0,
				0,
				0,
				0,
				0,
				0,
				0,
				0,
				-0.0333333333333333,
				0,
				0,
				0,
				0,
				0,
				0,
				0,
				0,
				0,
				0,
				0,
				0,
				0,
				0,
				0,
				0,
				0,
				0,
				0,
				0,
				0,
				-0.0333333333333333,
				0,
				0,
				0,
				0,
				0,
				0,
				0,
				0,
				0,
				0,
				0,
				0,
				0,
				0,
				0,
				0,
				0,
				0,
				0,
				0,
				0,
				-0.0333333333333333,
				0,
				0,
				0,
				0,
				0,
				0,
				0,
				0,
				0,
				0,
				0,
				0,
				0,
				0,
				0,
				0,
				0,
				0,
				0,
				0,
				0,
				-0.0333333333333333,
				0,
				0,
				0,
				0,
				0,
				0,
				0,
				0,
				0,
				0,
				0,
				0,
				0,
				0,
				0,
				0,
				0,
				0,
				0,
				0,
				0,
				-0.0333333333333333,
				0,
				0,
				0,
				0,
				0,
				0,
				0,
				0,
				0,
				0,
				0,
				0,
				0,
				0,
				0,
				0,
				0,
				0,
				0,
				0,
				0,
				-0.0250000000000000,
				0,
				0,
				0,
				0,
				0,
				0,
				0,
				0,
				0,
				0,
				0,
				0,
				0,
				0,
				0,
				0,
				0,
				0,
				0,
				0,
				0,
				-0.0166666666666667,
				-0.0166666666666667,
				0,
				-0.0250000000000000,
				0,
				-0.0333333333333333,
				0,
				-0.0333333333333333,
				0,
				-0.0333333333333333,
				0,
				-0.0333333333333333,
				0,
				-0.0333333333333333,
				0,
				-0.0333333333333333,
				0,
				-0.0250000000000000,
				0,
				-0.0166666666666667,
				0,
				-0.00833333333333334,
				-0.00833333333333334,
			};

			for (int i = 0; i < forceVectorExpected.Length; i++)
				Assert.Equal(forceVectorExpected[i], model.PatchesDictionary[0].Forces[i], 7);

			#region expectedDisplacement

			var displacementVectorExpected = new double[]
			{
				-0.00314820363909275,
				-0.00314820363909275,
				-0.00379395242120674,
				-0.00474046877963962,
				-0.00527896515562401,
				-0.00627080476725420,
				-0.00611202232103848,
				-0.00732618402014196,
				-0.00671664390713507,
				-0.00830517563434155,
				-0.00750778852140214,
				-0.00921569343748655,
				-0.00711656646035208,
				-0.0104552048087267,
				-0.00857783159759678,
				-0.0115132239445180,
				-0.00447714927679477,
				-0.0142326173251668,
				-0.00896507378072984,
				-0.0172326398522936,
				0.0165696344045362,
				-0.0344920417612453,
				-0.00474046877963962,
				-0.00379395242120674,
				-0.0106163667467970,
				-0.0106163667467970,
				-0.0145901385325659,
				-0.0159579381873314,
				-0.0175798662746167,
				-0.0200752820397931,
				-0.0197642416239734,
				-0.0232933594058257,
				-0.0209594659687415,
				-0.0269739141344681,
				-0.0221701240269120,
				-0.0304381009757610,
				-0.0204015380221106,
				-0.0372162332314755,
				-0.0190262809325074,
				-0.0451261112082839,
				0.000524263735174436,
				-0.0646122725733136,
				0.0134058371171140,
				-0.0674237539857168,
				-0.00627080476725419,
				-0.00527896515562400,
				-0.0159579381873313,
				-0.0145901385325660,
				-0.0255583394854834,
				-0.0255583394854834,
				-0.0320458295432382,
				-0.0340345779839745,
				-0.0364394696445877,
				-0.0416128738023956,
				-0.0392520264843806,
				-0.0489264479827333,
				-0.0394128785675202,
				-0.0579962012559210,
				-0.0375262954335112,
				-0.0685615165413566,
				-0.0269408733978005,
				-0.0856583155480079,
				-0.0140932604372300,
				-0.0958837290337786,
				-0.00425322389965576,
				-0.102864805751424,
				-0.00732618402014194,
				-0.00611202232103848,
				-0.0200752820397931,
				-0.0175798662746167,
				-0.0340345779839744,
				-0.0320458295432382,
				-0.0442399103415714,
				-0.0442399103415714,
				-0.0511182027891641,
				-0.0551752972826913,
				-0.0547664019599041,
				-0.0662589972721233,
				-0.0554544754480729,
				-0.0781539846420991,
				-0.0514092081933032,
				-0.0928498000772380,
				-0.0427832850822321,
				-0.108645057746992,
				-0.0311333714176029,
				-0.120963904936585,
				-0.0251089003632628,
				-0.125275523366600,
				-0.00830517563434156,
				-0.00671664390713508,
				-0.0232933594058256,
				-0.0197642416239734,
				-0.0416128738023955,
				-0.0364394696445877,
				-0.0551752972826911,
				-0.0511182027891642,
				-0.0644659012290671,
				-0.0644659012290672,
				-0.0694887264699124,
				-0.0777115338314980,
				-0.0702085098844087,
				-0.0920154063898269,
				-0.0670415920146972,
				-0.107342294443253,
				-0.0590775341097851,
				-0.124438661501351,
				-0.0522405165453920,
				-0.134720354624928,
				-0.0486371118189779,
				-0.139962781550140,
				-0.00921569343748652,
				-0.00750778852140214,
				-0.0269739141344680,
				-0.0209594659687415,
				-0.0489264479827331,
				-0.0392520264843806,
				-0.0662589972721231,
				-0.0547664019599042,
				-0.0777115338314978,
				-0.0694887264699125,
				-0.0840291771115493,
				-0.0840291771115496,
				-0.0855834304606439,
				-0.0992380743950709,
				-0.0829894406989148,
				-0.115516212797889,
				-0.0776932540445211,
				-0.131918682864287,
				-0.0732614094641817,
				-0.142841837190617,
				-0.0714730797892666,
				-0.147708704990261,
				-0.0104552048087266,
				-0.00711656646035210,
				-0.0304381009757609,
				-0.0221701240269120,
				-0.0579962012559208,
				-0.0394128785675203,
				-0.0781539846420988,
				-0.0554544754480730,
				-0.0920154063898266,
				-0.0702085098844090,
				-0.0992380743950705,
				-0.0855834304606441,
				-0.101501531630957,
				-0.101501531630958,
				-0.0999370644850315,
				-0.118017394845313,
				-0.0962343442201912,
				-0.134878482285728,
				-0.0940899630707035,
				-0.145633866412365,
				-0.0936977503670338,
				-0.150902312951907,
				-0.0115132239445179,
				-0.00857783159759677,
				-0.0372162332314754,
				-0.0204015380221106,
				-0.0685615165413564,
				-0.0375262954335112,
				-0.0928498000772377,
				-0.0514092081933034,
				-0.107342294443253,
				-0.0670415920146973,
				-0.115516212797889,
				-0.0829894406989151,
				-0.118017394845313,
				-0.0999370644850319,
				-0.117164741332239,
				-0.117164741332239,
				-0.114732961045519,
				-0.134270598345563,
				-0.113764345865824,
				-0.145437001617946,
				-0.114021716983238,
				-0.150889885103904,
				-0.0142326173251668,
				-0.00447714927679481,
				-0.0451261112082838,
				-0.0190262809325074,
				-0.0856583155480077,
				-0.0269408733978006,
				-0.108645057746992,
				-0.0427832850822322,
				-0.124438661501351,
				-0.0590775341097853,
				-0.131918682864286,
				-0.0776932540445212,
				-0.134878482285727,
				-0.0962343442201915,
				-0.134270598345562,
				-0.114732961045519,
				-0.132671391727535,
				-0.132671391727535,
				-0.132190010013491,
				-0.144083210808925,
				-0.132499710574820,
				-0.149832975007404,
				-0.0172326398522934,
				-0.00896507378072984,
				-0.0646122725733136,
				0.000524263735174192,
				-0.0958837290337782,
				-0.0140932604372301,
				-0.120963904936585,
				-0.0311333714176032,
				-0.134720354624928,
				-0.0522405165453921,
				-0.142841837190616,
				-0.0732614094641819,
				-0.145633866412364,
				-0.0940899630707035,
				-0.145437001617945,
				-0.113764345865824,
				-0.144083210808926,
				-0.132190010013491,
				-0.143845673350073,
				-0.143845673350073,
				-0.143883483258064,
				-0.149631002634067,
				-0.0344920417612453,
				0.0165696344045361,
				-0.0674237539857166,
				0.0134058371171139,
				-0.102864805751424,
				-0.00425322389965601,
				-0.125275523366599,
				-0.0251089003632630,
				-0.139962781550140,
				-0.0486371118189782,
				-0.147708704990260,
				-0.0714730797892668,
				-0.150902312951907,
				-0.0936977503670341,
				-0.150889885103902,
				-0.114021716983238,
				-0.149832975007405,
				-0.132499710574819,
				-0.149631002634066,
				-0.143883483258064,
				-0.149670908292985,
				-0.149670908292985,
			};

			for (int i = 0; i < displacementVectorExpected.Length; i++)
				Assert.Equal(displacementVectorExpected[i], solver.LinearSystems[0].Solution[i], 7);

			#endregion
		}

		//[Fact]
		public void IsogeometricPlateWithHole()
		{
			// Model
			VectorExtensions.AssignTotalAffinityCount();
			Model model = new Model();
			ModelCreator modelCreator = new ModelCreator(model);
			string filename = "..\\..\\..\\IGA\\InputFiles\\PlateTension.txt";
			IsogeometricReader modelReader = new IsogeometricReader(modelCreator, filename);
			modelReader.CreateModelFromFile();

			// Forces and Boundary Conditions
			Value horizontalDistributedLoad = delegate (double x, double y, double z)
			{
				return new double[] { 100, 0, 0 };
			};
			model.PatchesDictionary[0].EdgesDictionary[1].LoadingConditions
				.Add(new NeumannBoundaryCondition(horizontalDistributedLoad));

			// Boundary Conditions - Dirichlet
			foreach (ControlPoint controlPoint in model.PatchesDictionary[0].EdgesDictionary[0].ControlPointsDictionary
				.Values)
			{
				model.ControlPointsDictionary[controlPoint.ID].Constrains.Add(new Constraint() { DOF = DOFType.Y });
			}

			model.ControlPointsDictionary[0].Constrains.Add(new Constraint() {DOF = DOFType.Y});
			
			// Solvers
			var solverBuilder = new SuiteSparseSolver.Builder();
			solverBuilder.DofOrderer = new DofOrderer(
				new NodeMajorDofOrderingStrategy(), new NullReordering());
			ISolver_v2 solver = solverBuilder.BuildSolver(model);

			// Structural problem provider
			var provider = new ProblemStructural_v2(model, solver);

			// Linear static analysis
			var childAnalyzer = new LinearAnalyzer_v2(solver);
			var parentAnalyzer = new StaticAnalyzer_v2(model, solver, provider, childAnalyzer);

			// Run the analysis
			parentAnalyzer.Initialize();
			parentAnalyzer.Solve();

			double[] forceVectorExpected =
			{
			};
			for (int i = 0; i < forceVectorExpected.Length; i++)
				Assert.Equal(forceVectorExpected[i], model.PatchesDictionary[0].Forces[i], 8);

			#region expectedDisplacement

			double[] displacementVectorExpected =
			{
			};
			for (int i = 0; i < displacementVectorExpected.Length; i++)
				Assert.Equal(displacementVectorExpected[i], solver.LinearSystems[0].Solution[i], 4);

			#endregion
		}

		//[Fact]
		public void IsogeometricCurvedBeamBenchmark()
		{
			// Model
			VectorExtensions.AssignTotalAffinityCount();
			Model model = new Model();
			ModelCreator modelCreator = new ModelCreator(model);
			string filename = "..\\..\\..\\InputFiles\\CurvedBeam.txt";
			IsogeometricReader modelReader = new IsogeometricReader(modelCreator, filename);
			modelReader.CreateModelFromFile();

			// Forces and Boundary Conditions
			model.PatchesDictionary[0].EdgesDictionary[0].LoadingConditions
				.Add(new PressureBoundaryCondition(30000));

			// Boundary Conditions - Dirichlet
			foreach (ControlPoint controlPoint in model.PatchesDictionary[0].EdgesDictionary[0].ControlPointsDictionary
				.Values)
			{
				model.ControlPointsDictionary[controlPoint.ID].Constrains.Add(new Constraint() { DOF = DOFType.Y });
			}

			foreach (ControlPoint controlPoint in model.PatchesDictionary[0].EdgesDictionary[1].ControlPointsDictionary
				.Values)
			{
				model.ControlPointsDictionary[controlPoint.ID].Constrains.Add(new Constraint() { DOF = DOFType.X });
			}

			// Solvers
			var solverBuilder = new SuiteSparseSolver.Builder();
			solverBuilder.DofOrderer = new DofOrderer(
				new NodeMajorDofOrderingStrategy(), new NullReordering());
			ISolver_v2 solver = solverBuilder.BuildSolver(model);

			// Structural problem provider
			var provider = new ProblemStructural_v2(model, solver);

			// Linear static analysis
			var childAnalyzer = new LinearAnalyzer_v2(solver);
			var parentAnalyzer = new StaticAnalyzer_v2(model, solver, provider, childAnalyzer);

			// Run the analysis
			parentAnalyzer.Initialize();
			parentAnalyzer.Solve();

			Matrix<double> forceVectorExpected= MatlabReader.Read<double>("CurvedBeam.mat", "forceVector");
			
			for (int i = 0; i < forceVectorExpected.RowCount; i++)
				Assert.Equal(forceVectorExpected.At(i,0), model.PatchesDictionary[0].Forces[i], 7);

			Matrix<double> displacementVectorExpected = MatlabReader.Read<double>("CurvedBeam.mat", "forceVector");

			for (int i = 0; i < displacementVectorExpected.RowCount; i++)
				Assert.Equal(displacementVectorExpected.At(i,0), solver.LinearSystems[0].Solution[i], 7);
		}

		
		//[Fact]
		public void IsogeometricBeam3D()
		{
			VectorExtensions.AssignTotalAffinityCount();
			Model model = new Model();
			ModelCreator modelCreator = new ModelCreator(model);
			string filename = "Beam3D";
			string filepath= $"..\\..\\..\\InputFiles\\{filename}.txt";
			IsogeometricReader modelReader = new IsogeometricReader(modelCreator, filepath);
			modelReader.CreateModelFromFile();

			// Forces and Boundary Conditions
			foreach (ControlPoint controlPoint in model.PatchesDictionary[0].FacesDictionary[1].ControlPointsDictionary
				.Values)
			{
				model.Loads.Add(new Load()
				{
					Amount = -100,
					ControlPoint = model.ControlPoints[controlPoint.ID],
					DOF = DOFType.Z
				});
			}


			// Boundary Conditions - Dirichlet
			foreach (ControlPoint controlPoint in model.PatchesDictionary[0].FacesDictionary[0].ControlPointsDictionary
				.Values)
			{
				model.ControlPointsDictionary[controlPoint.ID].Constrains.Add(new Constraint() {DOF = DOFType.X});
				model.ControlPointsDictionary[controlPoint.ID].Constrains.Add(new Constraint() { DOF = DOFType.Y });
				model.ControlPointsDictionary[controlPoint.ID].Constrains.Add(new Constraint() { DOF = DOFType.Z });
			}

			// Solvers
			var solverBuilder = new SuiteSparseSolver.Builder();
			solverBuilder.DofOrderer = new DofOrderer(
				new NodeMajorDofOrderingStrategy(), new NullReordering());
			ISolver_v2 solver = solverBuilder.BuildSolver(model);

			// Structural problem provider
			var provider = new ProblemStructural_v2(model, solver);

			// Linear static analysis
			var childAnalyzer = new LinearAnalyzer_v2(solver);
			var parentAnalyzer = new StaticAnalyzer_v2(model, solver, provider, childAnalyzer);

			// Run the analysis
			parentAnalyzer.Initialize();
			parentAnalyzer.Solve();

			var paraview= new ParaviewNurbs3D(model, solver.LinearSystems[0].Solution, filename);
			paraview.CreateParaviewFile();
			//Matrix<double> forceVectorExpected = MatlabReader.Read<double>("..\\..\\..\\InputFiles\\Beam3D.mat", "forceVector");

			//for (int i = 0; i < forceVectorExpected.RowCount; i++)
			//	Assert.True(Utilities.AreValuesEqual(forceVectorExpected.At(i, 0), model.PatchesDictionary[0].Forces[i], 1e-2));

			//Matrix<double> displacementVectorExpected = MatlabReader.Read<double>("..\\..\\..\\InputFiles\\Beam3D.mat", "displacementVector");

			//for (int i = 0; i < displacementVectorExpected.RowCount; i++)
			//	Assert.True(Utilities.AreValuesEqual(displacementVectorExpected.At(i, 0), linearSystems[0].Solution[i], 1e-2));
		}

		//[Fact]
		public void TSplinesShellsBenchmark()
		{
			VectorExtensions.AssignTotalAffinityCount();
			Model model = new Model();
			string filename = "..\\..\\..\\InputFiles\\surface.iga";
			IGAFileReader modelReader = new IGAFileReader(model, filename);
			modelReader.CreateTSplineShellsModelFromFile();

			model.PatchesDictionary[0].Material = new ElasticMaterial2D(StressState2D.PlaneStress)
			{
				PoissonRatio = 0.3,
				YoungModulus = 10e6
			};
			model.PatchesDictionary[0].Thickness = 0.1;

			for (int i = 0; i < 100; i++)
			{
				model.ControlPointsDictionary[i].Constrains.Add(new Constraint() { DOF = DOFType.X });
				model.ControlPointsDictionary[i].Constrains.Add(new Constraint() { DOF = DOFType.Y });
				model.ControlPointsDictionary[i].Constrains.Add(new Constraint() { DOF = DOFType.Z });
			}

			for (int i = 0; i < model.ControlPoints.Count - 100; i++)
			{
				model.Loads.Add(new Load()
				{
					Amount = -1000,
					ControlPoint = model.ControlPointsDictionary[i],
					DOF = DOFType.Y
				});
			}


			// Solvers
			var solverBuilder = new SuiteSparseSolver.Builder();
			solverBuilder.DofOrderer = new DofOrderer(
				new NodeMajorDofOrderingStrategy(), new NullReordering());
			ISolver_v2 solver = solverBuilder.BuildSolver(model);

			// Structural problem provider
			var provider = new ProblemStructural_v2(model, solver);

			// Linear static analysis
			var childAnalyzer = new LinearAnalyzer_v2(solver);
			var parentAnalyzer = new StaticAnalyzer_v2(model, solver, provider, childAnalyzer);

			// Run the analysis
			parentAnalyzer.Initialize();
			parentAnalyzer.Solve();
		}
	}
}