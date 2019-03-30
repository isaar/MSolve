using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using ISAAR.MSolve.Discretization;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.FEM.Elements;
using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.Materials;
using ISAAR.MSolve.Preprocessor.Meshes;
using ISAAR.MSolve.Preprocessor.Meshes.Custom;

// Geometry:
//             | 
//             |
//            \|/ load
// |           . 
// |------------
// |    length

// Cross section:
//     |
//     | 
//    \|/ load
//     .
//   -----
//  |     |
//  |     |
//  |     | <--- height
//  |     |
//  |     |
//   -----
//     ^
//     |     
//   width

namespace ISAAR.MSolve.Solvers.Tests
{
    public class CantileverBeam
    {
        private const double eulerBeamDimensionsRatio = 10.0;
        private const int subdomainID = 0;

        private readonly double length;
        private readonly double height;
        private readonly double width;
        private readonly double endPointLoad;
        private readonly double youngModulus;
        private readonly Node_v2[] endNodes;

        private CantileverBeam(double length, double height, double width, double endPointLoad, 
            double youngModulus, Model_v2 model, Node_v2[] endNodes)
        {
            this.length = length;
            this.height = height;
            this.width = width;
            this.endPointLoad = endPointLoad;
            this.youngModulus = youngModulus;
            this.Model = model;
            this.endNodes = endNodes;
        }

        public Model_v2 Model { get; }

        public double CalculateEndDeflectionWithEulerBeamTheory()
        {
            if ((length < eulerBeamDimensionsRatio * height) || (length < eulerBeamDimensionsRatio * width))
            {
                throw new ArgumentException("In order for Euler beam theory to hold, the beam's length must be at least"
                    + $" {eulerBeamDimensionsRatio} times greater than each of the other two dimensions.");
            }
            double momentOfInertia = 1.0 / 12.0 * width * Math.Pow(height, 3.0);
            return endPointLoad * Math.Pow(length, 3.0) / (3.0 * youngModulus * momentOfInertia);
        }

        public double CalculateAverageEndDeflectionFromSolution(IVectorView solution)
        {
            //TODO: better do this with observers/loggers
            DofTable subdomainDofs = Model.SubdomainsDictionary[subdomainID].FreeDofOrdering.FreeDofs;
            double endDeflectionSum = 0.0;
            int dofsCount = 0;
            foreach (var node in endNodes)
            {
                bool exists = subdomainDofs.TryGetValue(node, DOFType.Y, out int dofIdx);
                if (exists)
                {
                    ++dofsCount;
                    endDeflectionSum += solution[dofIdx];
                }
            }
            return endDeflectionSum / dofsCount;
        }


        public class Builder
        {
            public double Length { get; set; } = 10.00;
            public double Height { get; set; } = 0.50;
            public double Width { get; set; } = 0.25;
            public double EndPointLoad { get; set; } = 20.0E3;
            public double YoungModulus { get; set; } = 2.1E7;
            public double PoissonRatio { get; set; } = 0.3;

            public CantileverBeam BuildWithBeamElements(int numElements)
            {
                throw new NotImplementedException();
            }

            public CantileverBeam BuildWithQuad4Elements(int numElementsAlongLength, int numElementsAlongHeight)
            {
                // Material and section properties
                double thickness = Width;
                var material = new ElasticMaterial2D_v2(StressState2D.PlaneStress)
                {
                    YoungModulus = this.YoungModulus,
                    PoissonRatio = this.PoissonRatio
                };
                var dynamicProperties = new DynamicMaterial(1.0, 0.0, 0.0);

                // Model with 1 subdomain
                var model = new Model_v2();
                model.SubdomainsDictionary.Add(subdomainID, new Subdomain_v2(subdomainID));

                // Generate mesh
                var meshGenerator = new UniformMeshGenerator2D_v2(0.0, 0.0, Length, Height, 
                    numElementsAlongLength, numElementsAlongHeight);
                (IReadOnlyList<Node_v2> vertices, IReadOnlyList<CellConnectivity_v2> cells) = meshGenerator.CreateMesh();

                // Add nodes to the model
                for (int n = 0; n < vertices.Count; ++n) model.NodesDictionary.Add(n, vertices[n]);

                // Add Quad4 elements to the model
                var factory = new ContinuumElement2DFactory(thickness, material, dynamicProperties);
                for (int e = 0; e < cells.Count; ++e)
                {
                    ContinuumElement2D element = factory.CreateElement(cells[e].CellType, cells[e].Vertices);
                    var elementWrapper = new Element_v2() { ID = e, ElementType = element };
                    foreach (Node_v2 node in element.Nodes) elementWrapper.AddNode(node);
                    model.ElementsDictionary.Add(e, elementWrapper);
                    model.SubdomainsDictionary[subdomainID].Elements.Add(elementWrapper);
                }

                // Clamp boundary condition at one end
                double tol = 1E-10; //TODO: this should be chosen w.r.t. the element size along X
                foreach (var node in model.Nodes.Where(node => Math.Abs(node.X) <= tol))
                {
                    node.Constraints.Add(new Constraint() { DOF = DOFType.X, Amount = 0.0 });
                    node.Constraints.Add(new Constraint() { DOF = DOFType.Y, Amount = 0.0 });
                }

                // Apply concentrated load at the other end
                Node_v2[] loadedNodes = model.Nodes.Where(node => Math.Abs(node.X - Length) <= tol).ToArray();
                foreach (var node in loadedNodes)
                {
                    model.Loads.Add(new Load_v2() { Amount = EndPointLoad / loadedNodes.Length, Node = node, DOF = DOFType.Y });
                }

                return new CantileverBeam(Length, Height, Width, EndPointLoad, YoungModulus, model, loadedNodes);
            }

            public CantileverBeam BuildWithHexa8Elements(int numElementsAlongLength, int numElementsAlongHeight, 
                int numElementsAlongWidth)
            {
                throw new NotImplementedException();
            }
        }
    }
}
