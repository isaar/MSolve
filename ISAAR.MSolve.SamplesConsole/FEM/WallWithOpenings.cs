using System;
using System.Collections.Generic;
using System.Linq;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.Discretization.Mesh.Generation;
using ISAAR.MSolve.Discretization.Mesh.Generation.GMSH;
using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.Materials;
using ISAAR.MSolve.SamplesConsole.Preprocessing;

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
            IReadOnlyList<Node> nodes;
            IReadOnlyList<CellConnectivity<Node>> elements;
            using (var reader = new GmshReader<Node>(meshPath))
            {
                (nodes, elements) = reader.CreateMesh((id, x, y, z) => new Node() { ID = id, X = x, Y = y, Z = z });
            }
            model.AddMesh2D(nodes, elements, material);

            // Prescribed displacements
            double tol = 1E-10;
            IEnumerable<Node> constrainedNodes = nodes.Where(node => Math.Abs(node.Y) <= tol);
            model.ApplyPrescribedDisplacements(constrainedNodes, StructuralDof.TranslationX, 0.0);
            model.ApplyPrescribedDisplacements(constrainedNodes, StructuralDof.TranslationY, 0.0);

            // Loads
            Node[] loadedNodes = nodes.Where(
                node => (Math.Abs(node.Y - height) <= tol) && ((Math.Abs(node.X) <= tol))).ToArray();
            if (loadedNodes.Length != 1) throw new Exception("Only 1 node was expected at the top left corner");
            model.ApplyNodalLoad(loadedNodes[0], StructuralDof.TranslationX, horizontalLoad);

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
            IReadOnlyList<Node> nodes;
            IReadOnlyList<CellConnectivity<Node>> elements;
            using (var reader = new GmshReader<Node>(meshPath))
            {
                (nodes, elements) = reader.CreateMesh((id, x, y, z) => new Node() { ID = id, X = x, Y = y, Z = z });
            }
            model.AddMesh2D(nodes, elements, material, dynamicProperties);

            // Prescribed displacements
            double tol = 1E-10;
            IEnumerable<Node> constrainedNodes = nodes.Where(node => Math.Abs(node.Y) <= tol);
            model.ApplyPrescribedDisplacements(constrainedNodes, StructuralDof.TranslationX, 0.0);
            model.ApplyPrescribedDisplacements(constrainedNodes, StructuralDof.TranslationY, 0.0);

            // Loads
            string accelerogramPath = workingDirectory + "\\elcentro_NS.txt";
            Dictionary<IDofType, double> magnifications = new Dictionary<IDofType, double>
            {
                { StructuralDof.TranslationX, 1.0 }
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