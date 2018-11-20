using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.IO;
using System.Linq;
using System.Text;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.FEM.Elements;
using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.Geometry.Shapes;
using ISAAR.MSolve.Logging.VTK;
using ISAAR.MSolve.Materials;
using ISAAR.MSolve.Preprocessor.Meshes;
using ISAAR.MSolve.Preprocessor.Meshes.Custom;
using ISAAR.MSolve.Preprocessor.Meshes.GMSH;
using ISAAR.MSolve.Preprocessor.UI;

namespace ISAAR.MSolve.SamplesConsole.FEM
{
    public class WallWithOpenings
    {
        public static void Run()
        {
            SolveStaticLinearWall();
            //SolveDynamicLinearWall();
        }

        private static void SolveStaticLinearWall()
        {
            // Some values
            string workingDirectory = @"C:\Users\Serafeim\Desktop\Presentation";
            double height = 3.5;
            double thickness = 0.1;
            double youngModulus = 2E6;
            double poissonRatio = 0.3;
            double horizontalLoad = 1000.0;

            // Initialize model
            PreprocessorModel model = PreprocessorModel.Create2DPlaneStress(thickness);

            // Materials
            ElasticMaterial2D material = new ElasticMaterial2D(StressState2D.PlaneStress)
            {
                YoungModulus = youngModulus,
                PoissonRatio = poissonRatio
            };

            // Read mesh from GMSH file
            string meshPath = workingDirectory + "\\wall.msh";
            IReadOnlyList<Node2D> nodes;
            IReadOnlyList<CellConnectivity2D> elements;
            using (var reader = new GmshReader2D(meshPath))
            {
                (nodes, elements) = reader.CreateMesh();
            }
            model.AddMesh2D(nodes, elements, material);

            // Prescribed displacements
            double tol = 1E-10;
            IEnumerable<Node2D> constrainedNodes = nodes.Where(node => Math.Abs(node.Y) <= tol);
            model.ApplyPrescribedDisplacements(constrainedNodes, DOFType.X, 0.0);
            model.ApplyPrescribedDisplacements(constrainedNodes, DOFType.Y, 0.0);

            // Loads
            Node2D[] loadedNodes = nodes.Where(
                node => (Math.Abs(node.Y - height) <= tol) && ((Math.Abs(node.X) <= tol))).ToArray();
            if (loadedNodes.Length != 1) throw new Exception("Only 1 node was expected at the top left corner");
            model.ApplyNodalLoad(loadedNodes[0], DOFType.X, horizontalLoad);

            // Define output
            OutputRequests output = new OutputRequests(workingDirectory + "\\Plots");
            output.Displacements = true;
            output.Strains = true;
            output.Stresses = true;
            output.StressesVonMises = true;

            // Set up the simulation procedure
            Job job = new Job(model);
            job.Procedure = Job.ProcedureOptions.Static;
            job.Integrator = Job.IntegratorOptions.Linear;
            job.Solver = Job.SolverOptions.DirectSkyline;
            job.FieldOutputRequests = output;

            // Run the simulation
            job.Submit();
        }

        private static void SolveDynamicLinearWall()
        {
            // Some values
            string workingDirectory = @"C:\Users\Serafeim\Desktop\Presentation";
            double thickness = 0.1;
            double youngModulus = 2E6;
            double poissonRatio = 0.3;

            // Initialize model
            PreprocessorModel model = PreprocessorModel.Create2DPlaneStress(thickness);

            // Materials
            ElasticMaterial2D material = new ElasticMaterial2D(StressState2D.PlaneStress)
            {
                YoungModulus = youngModulus,
                PoissonRatio = poissonRatio
            };
            DynamicMaterial dynamicProperties = new DynamicMaterial(25, 0.05, 0.05);

            // Read mesh from GMSH file
            string meshPath = workingDirectory + "\\wall.msh";
            IReadOnlyList<Node2D> nodes;
            IReadOnlyList<CellConnectivity2D> elements;
            using (var reader = new GmshReader2D(meshPath))
            {
                (nodes, elements) = reader.CreateMesh();
            }
            model.AddMesh2D(nodes, elements, material, dynamicProperties);

            // Prescribed displacements
            double tol = 1E-10;
            IEnumerable<Node2D> constrainedNodes = nodes.Where(node => Math.Abs(node.Y) <= tol);
            model.ApplyPrescribedDisplacements(constrainedNodes, DOFType.X, 0.0);
            model.ApplyPrescribedDisplacements(constrainedNodes, DOFType.Y, 0.0);

            // Loads
            string accelerogramPath = workingDirectory + "\\elcentro_NS.txt";
            Dictionary<DOFType, double> magnifications = new Dictionary<DOFType, double>
            {
                { DOFType.X, 1.0 }
            };
            model.SetGroundMotion(accelerogramPath, magnifications, 0.02, 53.74);

            // Define output
            OutputRequests output = new OutputRequests(workingDirectory + "\\Plots");
            output.Displacements = true;
            output.Strains = true;
            output.Stresses = true;
            output.StressesVonMises = true;

            // Set up the simulation procedure
            Job job = new Job(model);
            job.Procedure = Job.ProcedureOptions.DynamicImplicit;
            job.Integrator = Job.IntegratorOptions.Linear;
            job.Solver = Job.SolverOptions.DirectSkyline;
            job.FieldOutputRequests = output;

            // Run the simulation
            job.Submit();
        }
    }
}