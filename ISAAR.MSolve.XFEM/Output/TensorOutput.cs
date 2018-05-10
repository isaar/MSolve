using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.XFEM.Elements;
using ISAAR.MSolve.XFEM.Entities;
using ISAAR.MSolve.XFEM.FreedomDegrees.Ordering;
using ISAAR.MSolve.XFEM.Geometry.CoordinateSystems;
using ISAAR.MSolve.XFEM.Interpolation.GaussPointSystems;
using ISAAR.MSolve.XFEM.Tensors;

namespace ISAAR.MSolve.XFEM.Output
{
    // TODO: This class should be more flexible, allowing the caller to specify which output he wants, but still not 
    // repeating calculations. The command pattern seems suitable to handle what fields will be calculated, 
    // extrapolated, etc and how
    class TensorOutput
    {
        private readonly Model2D model;
        private readonly IDofOrderer dofOrderer;

        public TensorOutput(Model2D model, IDofOrderer dofOrderer)
        {
            this.model = model;
            this.dofOrderer = dofOrderer;
        }

        public IReadOnlyDictionary<XContinuumElement2D, IReadOnlyList<Tensor2D>> ComputeElementWiseNodalStrains(
            Vector solution, bool extrapolateFromGPs = false)
        {
            return ComputeElementWiseNodalTensors(solution, extrapolateFromGPs, new StrainField());
        }

        public IReadOnlyDictionary<XContinuumElement2D, IReadOnlyList<Tensor2D>> ComputeElementWiseNodalStresses(
            Vector solution, bool extrapolateFromGPs = false)
        {
            return ComputeElementWiseNodalTensors(solution, extrapolateFromGPs, new StressField());
        }

        public IReadOnlyList<Tensor2D> ComputeSmoothedNodalStrains(Vector solution, bool extrapolateFromGPs = false)
        {
            return ComputeSmoothedNodalTensors(solution, extrapolateFromGPs, new StrainField());
        }

        public IReadOnlyList<Tensor2D> ComputeSmoothedNodalStresses(Vector solution, bool extrapolateFromGPs = false)
        {
            return ComputeSmoothedNodalTensors(solution, extrapolateFromGPs, new StressField());
        }

        private IReadOnlyList<Tensor2D> ComputeSmoothedNodalTensors(Vector solution, bool extrapolateFromGPs, 
            IOutputField field)
        {
            var tensorsFromAllElements = new Dictionary<XNode2D, List<Tensor2D>>();
            foreach (var node in model.Nodes) tensorsFromAllElements[node] = new List<Tensor2D>();
            Vector constrainedDisplacements = model.CalculateConstrainedDisplacements(dofOrderer);

            foreach (var element in model.Elements)
            {
                IReadOnlyList<Tensor2D> elementTensors;
                if (extrapolateFromGPs)
                {
                    elementTensors = 
                        ElementNodalTensorsExtrapolation(element, solution, constrainedDisplacements, field);
                }
                else
                {
                    elementTensors = ElementNodalTensorsDirectly(element, solution, constrainedDisplacements, field);
                }
                
                for (int nodeIdx = 0; nodeIdx < element.Nodes.Count; ++nodeIdx)
                {
                    tensorsFromAllElements[element.Nodes[nodeIdx]].Add(elementTensors[nodeIdx]);
                }
            }

            // Average with equal weights for all elements. TODO: perhaps vary the weights depending on the element type/area
            var nodalTensors = new Tensor2D[model.Nodes.Count];
            for (int i = 0; i < model.Nodes.Count; ++i)
            {
                XNode2D node = model.Nodes[i];
                double tensorXX = 0.0, tensorYY = 0.0, tensorXY = 0.0;
                foreach (var tensor in tensorsFromAllElements[node])
                {
                    tensorXX += tensor.XX;
                    tensorYY += tensor.YY;
                    tensorXY += tensor.XY;
                }
                int contributingElementsCount = tensorsFromAllElements[node].Count;
                tensorXX /= contributingElementsCount;
                tensorYY /= contributingElementsCount;
                tensorXY /= contributingElementsCount;
                nodalTensors[i] = new Tensor2D(tensorXX, tensorYY, tensorXY);
            }

            return nodalTensors;
        }

        private IReadOnlyDictionary<XContinuumElement2D, IReadOnlyList<Tensor2D>> ComputeElementWiseNodalTensors(
            Vector solution, bool extrapolateFromGPs, IOutputField field)
        {
            Vector constrainedDisplacements = model.CalculateConstrainedDisplacements(dofOrderer);
            var allTensors = new Dictionary<XContinuumElement2D, IReadOnlyList<Tensor2D>>();
            foreach (var element in model.Elements)
            {
                if (extrapolateFromGPs)
                {
                    allTensors[element] =
                        ElementNodalTensorsExtrapolation(element, solution, constrainedDisplacements, field);
                }
                else
                {
                    allTensors[element] =
                        ElementNodalTensorsDirectly(element, solution, constrainedDisplacements, field);
                }
            }
            return allTensors;
        }

        // Computes stresses directly at the nodes. The other approach is to compute them at Gauss points and then extrapolate
        private IReadOnlyList<Tensor2D> ElementNodalTensorsDirectly(XContinuumElement2D element,
            Vector freeDisplacements, Vector constrainedDisplacements, IOutputField field)
        {
            Vector standardDisplacements = dofOrderer.ExtractDisplacementVectorOfElementFromGlobal(element,
                freeDisplacements, constrainedDisplacements);
            Vector enrichedDisplacements =
                dofOrderer.ExtractEnrichedDisplacementsOfElementFromGlobal(element, freeDisplacements);

            IReadOnlyList<INaturalPoint2D> naturalNodes = element.ElementType.NaturalCoordinatesOfNodes;
            var nodalTensors = new Tensor2D[element.Nodes.Count];
            for (int i = 0; i < element.Nodes.Count; ++i)
            {
                nodalTensors[i] =
                    field.EvaluateAt(element, naturalNodes[i], standardDisplacements, enrichedDisplacements);
            }
            return nodalTensors;
        }

        private IReadOnlyList<Tensor2D> ElementNodalTensorsExtrapolation(XContinuumElement2D element,
            Vector freeDisplacements, Vector constrainedDisplacements, IOutputField field)
        {
            Vector standardDisplacements = dofOrderer.ExtractDisplacementVectorOfElementFromGlobal(element,
                freeDisplacements, constrainedDisplacements);
            Vector enrichedDisplacements =
                dofOrderer.ExtractEnrichedDisplacementsOfElementFromGlobal(element, freeDisplacements);

            IGaussPointSystem gpSystem = element.ElementType.GaussPointSystem;
            IReadOnlyList<INaturalPoint2D> gaussPoints = gpSystem.GaussPoints;
            var gpTensors = new Tensor2D[gaussPoints.Count];
            for (int i = 0; i < gaussPoints.Count; ++i)
            {
                gpTensors[i] =
                    field.EvaluateAt(element, gaussPoints[i], standardDisplacements, enrichedDisplacements);
            }
            return gpSystem.ExtrapolateTensorFromGaussPointsToNodes(gpTensors);
        }

    }
}
