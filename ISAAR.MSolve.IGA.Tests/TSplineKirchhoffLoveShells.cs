using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using ISAAR.MSolve.Analyzers;
using ISAAR.MSolve.Discretization;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.IGA.Entities;
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
using Xunit;

namespace ISAAR.MSolve.IGA.Tests
{
	public class TSplineKirchhoffLoveShells
	{
		[Fact]
		public void CantileverShellBenchmark()
		{
			VectorExtensions.AssignTotalAffinityCount();
			Model model = new Model();
			var filename = "CantileverShell";
			string filepath = $"..\\..\\..\\InputFiles\\{filename}.iga";
			IGAFileReader modelReader = new IGAFileReader(model, filepath);
			modelReader.CreateTSplineShellsModelFromFile();

			model.PatchesDictionary[0].Material = new ElasticMaterial2D(StressState2D.PlaneStress)
			{
				PoissonRatio = 0.0,
				YoungModulus = 100
			};
			model.PatchesDictionary[0].Thickness = 1;

			foreach (var controlPoint in model.ControlPointsDictionary.Values.Where(cp=>cp.X<3))
			{
				model.ControlPointsDictionary[controlPoint.ID].Constrains.Add(new Constraint(){ DOF = DOFType.X});
				model.ControlPointsDictionary[controlPoint.ID].Constrains.Add(new Constraint() { DOF = DOFType.Y });
				model.ControlPointsDictionary[controlPoint.ID].Constrains.Add(new Constraint() { DOF = DOFType.Z });
			}

			foreach (var controlPoint in model.ControlPointsDictionary.Values.Where(cp => cp.X >49.8))
			{
				model.Loads.Add(new Load()
				{
					Amount = -0.5,
					ControlPoint = model.ControlPointsDictionary[controlPoint.ID],
					DOF = DOFType.Z
				});
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

			var paraview = new ParaviewTsplineShells(model, solver.LinearSystems[0].Solution, filename);
			paraview.CreateParaviewFile();

			var expectedSolutionVector = new Vector(new double[]
			{
				0, 0, -306.122431, 0, 0, -1552.478121, 0, 0, -3454.810388, 0, 0, -5881.924153, 0, 0, -8702.62361, 0, 0,
				-11785.71439, 0, 0, -13928.57064, 0, 0, -15000.0008, 0, 0, -306.1224369, 0, 0, -1552.47811, 0, 0,
				-3454.810407, 0, 0, -5881.924117, 0, 0, -8702.623683, 0, 0, -11785.71423, 0, 0, -13928.57093, 0, 0,
				-15000.00025, 0, 0, -306.1224493, 0, 0, -1552.478088, 0, 0, -3454.810449, 0, 0, -5881.924038, 0, 0,
				-8702.623837, 0, 0, -11785.71389, 0, 0, -13928.57157, 0, 0, -14999.99909, 0, 0, -306.1224494, 0, 0,
				-1552.478088, 0, 0, -3454.810449, 0, 0, -5881.924038, 0, 0, -8702.623837, 0, 0, -11785.71389, 0, 0,
				-13928.57157, 0, 0, -14999.99909, 0, 0, -306.1224369, 0, 0, -1552.47811, 0, 0, -3454.810407, 0, 0,
				-5881.924117, 0, 0, -8702.623683, 0, 0, -11785.71423, 0, 0, -13928.57093, 0, 0, -15000.00025, 0, 0,
				-306.122431, 0, 0, -1552.478121, 0, 0, -3454.810388, 0, 0, -5881.924154, 0, 0, -8702.62361, 0, 0,
				-11785.71439, 0, 0, -13928.57064, 0, 0, -15000.0008
			});
			for (int i = 0; i < expectedSolutionVector.Length; i++)
			{
				Assert.True(Utilities.AreValuesEqual(expectedSolutionVector[i], solver.LinearSystems[0].Solution[i],
					1e-8));
			}
			
		}


		[Fact]
		public void CantileverShellMaterialBenchmark()
		{
			VectorExtensions.AssignTotalAffinityCount();
			Model model = new Model();
			string filename = "..\\..\\..\\InputFiles\\CantileverShell.iga";
			IGAFileReader modelReader = new IGAFileReader(model, filename);
			
			var thickness = 1.0;

			modelReader.CreateTSplineShellsModelFromFile(IGAFileReader.TSplineShellTypes.ThicknessMaterial, new ShellElasticMaterial2D
			{
				PoissonRatio = 0.0,
				YoungModulus = 100,
			}, thickness);
			foreach (var controlPoint in model.ControlPointsDictionary.Values.Where(cp => cp.X < 3))
			{
				model.ControlPointsDictionary[controlPoint.ID].Constrains.Add(new Constraint(){DOF = DOFType.X});
				model.ControlPointsDictionary[controlPoint.ID].Constrains.Add(new Constraint() { DOF = DOFType.Y });
				model.ControlPointsDictionary[controlPoint.ID].Constrains.Add(new Constraint() { DOF = DOFType.Z });
			}

			foreach (var controlPoint in model.ControlPointsDictionary.Values.Where(cp => cp.X > 49.8))
			{
				model.Loads.Add(new Load()
				{
					Amount = -0.5,
					ControlPoint = model.ControlPointsDictionary[controlPoint.ID],
					DOF = DOFType.Z
				});
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

			var expectedSolutionVector = new Vector(new double[]
			{
				0, 0, -306.122431, 0, 0, -1552.478121, 0, 0, -3454.810388, 0, 0, -5881.924153, 0, 0, -8702.62361, 0, 0,
				-11785.71439, 0, 0, -13928.57064, 0, 0, -15000.0008, 0, 0, -306.1224369, 0, 0, -1552.47811, 0, 0,
				-3454.810407, 0, 0, -5881.924117, 0, 0, -8702.623683, 0, 0, -11785.71423, 0, 0, -13928.57093, 0, 0,
				-15000.00025, 0, 0, -306.1224493, 0, 0, -1552.478088, 0, 0, -3454.810449, 0, 0, -5881.924038, 0, 0,
				-8702.623837, 0, 0, -11785.71389, 0, 0, -13928.57157, 0, 0, -14999.99909, 0, 0, -306.1224494, 0, 0,
				-1552.478088, 0, 0, -3454.810449, 0, 0, -5881.924038, 0, 0, -8702.623837, 0, 0, -11785.71389, 0, 0,
				-13928.57157, 0, 0, -14999.99909, 0, 0, -306.1224369, 0, 0, -1552.47811, 0, 0, -3454.810407, 0, 0,
				-5881.924117, 0, 0, -8702.623683, 0, 0, -11785.71423, 0, 0, -13928.57093, 0, 0, -15000.00025, 0, 0,
				-306.122431, 0, 0, -1552.478121, 0, 0, -3454.810388, 0, 0, -5881.924154, 0, 0, -8702.62361, 0, 0,
				-11785.71439, 0, 0, -13928.57064, 0, 0, -15000.0008
			});
			for (int i = 0; i < expectedSolutionVector.Length; i++)
			{
				Assert.True(Utilities.AreValuesEqual(expectedSolutionVector[i], solver.LinearSystems[0].Solution[i],
					1e-6));
			}

		}

		[Fact]
		public void CantileverShellMaterialBenchmark_v2()
		{
			VectorExtensions.AssignTotalAffinityCount();
			Model model = new Model();
			string filename = "..\\..\\..\\InputFiles\\CantileverShell.iga";
			IGAFileReader modelReader = new IGAFileReader(model, filename);

			var thickness = 1.0;

			modelReader.CreateTSplineShellsModelFromFile(IGAFileReader.TSplineShellTypes.ThicknessMaterial, new ShellElasticMaterial2D
			{
				PoissonRatio = 0.0,
				YoungModulus = 100,
			}, thickness);
			foreach (var controlPoint in model.ControlPointsDictionary.Values.Where(cp => cp.X < 3))
			{
				model.ControlPointsDictionary[controlPoint.ID].Constrains.Add(new Constraint(){DOF = DOFType.X});
				model.ControlPointsDictionary[controlPoint.ID].Constrains.Add(new Constraint() { DOF = DOFType.Y });
				model.ControlPointsDictionary[controlPoint.ID].Constrains.Add(new Constraint() { DOF = DOFType.Z });
			}

			foreach (var controlPoint in model.ControlPointsDictionary.Values.Where(cp => cp.X > 49.8))
			{
				model.Loads.Add(new Load()
				{
					Amount = -0.5,
					ControlPoint = model.ControlPointsDictionary[controlPoint.ID],
					DOF = DOFType.Z
				});
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


			var expectedSolutionVector = new Vector(new double[]
			{
				0, 0, -306.122431, 0, 0, -1552.478121, 0, 0, -3454.810388, 0, 0, -5881.924153, 0, 0, -8702.62361, 0, 0,
				-11785.71439, 0, 0, -13928.57064, 0, 0, -15000.0008, 0, 0, -306.1224369, 0, 0, -1552.47811, 0, 0,
				-3454.810407, 0, 0, -5881.924117, 0, 0, -8702.623683, 0, 0, -11785.71423, 0, 0, -13928.57093, 0, 0,
				-15000.00025, 0, 0, -306.1224493, 0, 0, -1552.478088, 0, 0, -3454.810449, 0, 0, -5881.924038, 0, 0,
				-8702.623837, 0, 0, -11785.71389, 0, 0, -13928.57157, 0, 0, -14999.99909, 0, 0, -306.1224494, 0, 0,
				-1552.478088, 0, 0, -3454.810449, 0, 0, -5881.924038, 0, 0, -8702.623837, 0, 0, -11785.71389, 0, 0,
				-13928.57157, 0, 0, -14999.99909, 0, 0, -306.1224369, 0, 0, -1552.47811, 0, 0, -3454.810407, 0, 0,
				-5881.924117, 0, 0, -8702.623683, 0, 0, -11785.71423, 0, 0, -13928.57093, 0, 0, -15000.00025, 0, 0,
				-306.122431, 0, 0, -1552.478121, 0, 0, -3454.810388, 0, 0, -5881.924154, 0, 0, -8702.62361, 0, 0,
				-11785.71439, 0, 0, -13928.57064, 0, 0, -15000.0008
			});
			for (int i = 0; i < expectedSolutionVector.Length; i++)
			{
				Assert.True(Utilities.AreValuesEqual(expectedSolutionVector[i], solver.LinearSystems[0].Solution[i],
					1e-6));
			}

		}

		//[Fact]
		public void SimpleHoodBenchmark()
		{
			VectorExtensions.AssignTotalAffinityCount();
			Model model = new Model();
			var filename = "attempt2";
			string filepath = $"..\\..\\..\\InputFiles\\{filename}.iga";
			IGAFileReader modelReader = new IGAFileReader(model, filepath);

			//var thickness = 1.0;

			//modelReader.CreateTSplineShellsModelFromFile(IGAFileReader.TSplineShellTypes.LinearMaterial,new ShellElasticMaterial2D
			//{
			//	PoissonRatio = 0.3,
			//	YoungModulus = 1e5,
			//}, thickness);
			modelReader.CreateTSplineShellsModelFromFile();

			model.PatchesDictionary[0].Material = new ElasticMaterial2D(StressState2D.PlaneStress)
			{
				PoissonRatio = 0.3,
				YoungModulus = 10000
			};
			model.PatchesDictionary[0].Thickness = 1;

			for (int i = 0; i < 100; i++)
			{
				var id = model.ControlPoints[i].ID;
				model.ControlPointsDictionary[id].Constrains.Add(new Constraint(){DOF = DOFType.X});
				model.ControlPointsDictionary[id].Constrains.Add(new Constraint() { DOF = DOFType.Y });
				model.ControlPointsDictionary[id].Constrains.Add(new Constraint() { DOF = DOFType.Z });
			}

			for (int i = model.ControlPoints.Count-100; i < model.ControlPoints.Count; i++)
			{
				var id = model.ControlPoints[i].ID;
				model.Loads.Add(new Load()
				{
					Amount = 100,
					ControlPoint = model.ControlPointsDictionary[id],
					DOF = DOFType.Z
				});
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

			var paraview= new ParaviewTsplineShells(model, solver.LinearSystems[0].Solution,filename);
			paraview.CreateParaviewFile();
		}

		[Fact]
		public void SimpleHoodBenchmarkMKL()
		{
			VectorExtensions.AssignTotalAffinityCount();
			Model model = new Model();
			var filename = "attempt2";
			string filepath = $"..\\..\\..\\InputFiles\\{filename}.iga";
			IGAFileReader modelReader = new IGAFileReader(model, filepath);

			var thickness = 1.0;
			
			modelReader.CreateTSplineShellsModelFromFile(IGAFileReader.TSplineShellTypes.ThicknessMaterial,new ShellElasticMaterial2D()
			{
				PoissonRatio = 0.3,
				YoungModulus = 10000
			}, thickness);
			

			for (int i = 0; i < 100; i++)
			{
				var id = model.ControlPoints[i].ID;
				model.ControlPointsDictionary[id].Constrains.Add(new Constraint(){DOF = DOFType.X});
				model.ControlPointsDictionary[id].Constrains.Add(new Constraint() { DOF = DOFType.Y });
				model.ControlPointsDictionary[id].Constrains.Add(new Constraint() { DOF = DOFType.Z });
			}

			for (int i = model.ControlPoints.Count - 100; i < model.ControlPoints.Count; i++)
			{
				var id = model.ControlPoints[i].ID;
				model.Loads.Add(new Load()
				{
					Amount = 100,
					ControlPoint = model.ControlPointsDictionary[id],
					DOF = DOFType.Z
				});
			}
			var solverBuilder = new SuiteSparseSolver.Builder();
			solverBuilder.DofOrderer = new DofOrderer(
				new NodeMajorDofOrderingStrategy(), AmdReordering.CreateWithSuiteSparseAmd());
			ISolver_v2 solver = solverBuilder.BuildSolver(model);

			// Structural problem provider
			var provider = new ProblemStructural_v2(model, solver);

			// Linear static analysis
			var childAnalyzer = new LinearAnalyzer_v2(solver);
			var parentAnalyzer = new StaticAnalyzer_v2(model, solver, provider, childAnalyzer);

			// Run the analysis
			parentAnalyzer.Initialize();
			parentAnalyzer.Solve();

			//var paraview = new ParaviewTsplineShells(model, solver.LinearSystems[0].Solution, filename);
			//paraview.CreateParaviewFile();
		}
	}
}
