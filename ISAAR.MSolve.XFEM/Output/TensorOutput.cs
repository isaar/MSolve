using System.Collections.Generic;
using System.Linq;
using ISAAR.MSolve.Discretization.Commons;
using ISAAR.MSolve.Geometry.Coordinates;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.XFEM.Elements;
using ISAAR.MSolve.XFEM.Entities;
using ISAAR.MSolve.XFEM.FreedomDegrees.Ordering;

namespace ISAAR.MSolve.XFEM.Output
{
    // TODO: This class should be more flexible, allowing the caller to specify which output he wants, but still not 
    // repeating calculations. The command pattern seems suitable to handle what fields will be calculated, 
    // extrapolated, etc and how
    class TensorOutput
    {
        private readonly Model2D_old model;
        private readonly IDofOrderer dofOrderer;

        public TensorOutput(Model2D_old model, IDofOrderer dofOrderer)
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
            var tensorsFromAllElements = new Dictionary<XNode, List<Tensor2D>>();
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
                XNode node = model.Nodes[i];
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

            IReadOnlyList<NaturalPoint> naturalNodes = element.Interpolation.NodalNaturalCoordinates;
            var nodalTensors = new Tensor2D[element.Nodes.Count];
            for (int i = 0; i < element.Nodes.Count; ++i)
            {
                nodalTensors[i] =
                    field.EvaluateAt(element, naturalNodes[i], standardDisplacements, enrichedDisplacements);
            }
            return nodalTensors;
        }

        //TODO: Either work with Tensor class or double[]
        private IReadOnlyList<Tensor2D> ElementNodalTensorsExtrapolation(XContinuumElement2D element,
            Vector freeDisplacements, Vector constrainedDisplacements, IOutputField field)
        {
            Vector standardDisplacements = dofOrderer.ExtractDisplacementVectorOfElementFromGlobal(element,
                freeDisplacements, constrainedDisplacements);
            Vector enrichedDisplacements =
                dofOrderer.ExtractEnrichedDisplacementsOfElementFromGlobal(element, freeDisplacements);

            IReadOnlyList<NaturalPoint> gaussPoints = element.GaussPointExtrapolation.Quadrature.IntegrationPoints;
            var gpTensors = new Tensor2D[gaussPoints.Count];
            for (int i = 0; i < gaussPoints.Count; ++i)
            {
                gpTensors[i] =
                    field.EvaluateAt(element, gaussPoints[i], standardDisplacements, enrichedDisplacements);
            }
            return element.GaussPointExtrapolation.ExtrapolateTensorFromGaussPointsToNodes(gpTensors, element.Interpolation);
        }

    }
}
