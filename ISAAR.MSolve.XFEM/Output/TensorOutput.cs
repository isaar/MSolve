using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.Numerical.LinearAlgebra;
using ISAAR.MSolve.XFEM.Elements;
using ISAAR.MSolve.XFEM.Entities;
using ISAAR.MSolve.XFEM.Geometry.CoordinateSystems;
using ISAAR.MSolve.XFEM.Interpolation;
using ISAAR.MSolve.XFEM.LinearAlgebra;
using ISAAR.MSolve.XFEM.Tensors;

namespace ISAAR.MSolve.XFEM.Output
{
    // TODO: This class should be more flexible, allowing the caller to specify which output he wants, but still not 
    // repeating calculations. The command pattern seems suitable to handle what fields will be calculated, 
    // extrapolated, etc and how
    class TensorOutput
    {
        private readonly Model2D model;

        public TensorOutput(Model2D model)
        {
            this.model = model;
        }

        public IReadOnlyDictionary<XContinuumElement2D, IReadOnlyList<Tensor2D>> ComputeElementWiseNodalStrains(
            double[] solution, bool extrapolateFromGPs = false)
        {
            return ComputeElementWiseNodalTensors(solution, extrapolateFromGPs, new StrainField());
        }

        public IReadOnlyDictionary<XContinuumElement2D, IReadOnlyList<Tensor2D>> ComputeElementWiseNodalStresses(
            double[] solution, bool extrapolateFromGPs = false)
        {
            return ComputeElementWiseNodalTensors(solution, extrapolateFromGPs, new StressField());
        }

        public IReadOnlyList<Tensor2D> ComputeSmoothedNodalStrains(double[] solution, bool extrapolateFromGPs = false)
        {
            return ComputeSmoothedNodalTensors(solution, extrapolateFromGPs, new StrainField());
        }

        public IReadOnlyList<Tensor2D> ComputeSmoothedNodalStresses(double[] solution, bool extrapolateFromGPs = false)
        {
            return ComputeSmoothedNodalTensors(solution, extrapolateFromGPs, new StressField());
        }

        private IReadOnlyList<Tensor2D> ComputeSmoothedNodalTensors(double[] solution, bool extrapolateFromGPs, 
            IOutputField field)
        {
            var tensorsFromAllElements = new Dictionary<XNode2D, List<Tensor2D>>();
            foreach (var node in model.Nodes) tensorsFromAllElements[node] = new List<Tensor2D>();
            double[] constrainedDisplacements = model.CalculateConstrainedDisplacements();

            foreach (var element in model.Elements)
            {
                IReadOnlyDictionary<XNode2D, Tensor2D> elementTensors;
                if (extrapolateFromGPs)
                {
                    elementTensors = 
                        ElementNodalTensorsExtrapolation(element, solution, constrainedDisplacements, field);
                }
                else
                {
                    elementTensors = ElementNodalTensorsDirectly(element, solution, constrainedDisplacements, field);
                }
                
                foreach (var nodeTensorPair in elementTensors)
                {
                    tensorsFromAllElements[nodeTensorPair.Key].Add(nodeTensorPair.Value);
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
            double[] solution, bool extrapolateFromGPs, IOutputField field)
        {
            double[] constrainedDisplacements = model.CalculateConstrainedDisplacements();
            var allTensors = new Dictionary<XContinuumElement2D, IReadOnlyList<Tensor2D>>();
            foreach (var element in model.Elements)
            {
                IReadOnlyDictionary<XNode2D, Tensor2D> elementTensors;
                if (extrapolateFromGPs)
                {
                    elementTensors = 
                        ElementNodalTensorsExtrapolation(element, solution, constrainedDisplacements, field);
                }
                else
                {
                    elementTensors = ElementNodalTensorsDirectly(element, solution, constrainedDisplacements, field);
                }
                allTensors[element] = elementTensors.Values.ToArray();
            }
            return allTensors;
        }

        // Computes stresses directly at the nodes. The other approach is to compute them at Gauss points and then extrapolate
        private IReadOnlyDictionary<XNode2D, Tensor2D> ElementNodalTensorsDirectly(XContinuumElement2D element,
            double[] freeDisplacements, double[] constrainedDisplacements, IOutputField field)
        {
            double[] standardDisplacements = model.DofEnumerator.ExtractDisplacementVectorOfElementFromGlobal(element,
                freeDisplacements, constrainedDisplacements);
            double[] enrichedDisplacements =
                model.DofEnumerator.ExtractEnrichedDisplacementsOfElementFromGlobal(element, freeDisplacements);

            IReadOnlyList<INaturalPoint2D> naturalNodes = element.ElementType.NaturalCoordinatesOfNodes;
            var nodalTensors = new Dictionary<XNode2D, Tensor2D>();
            for (int i = 0; i < element.Nodes.Count; ++i)
            {
                nodalTensors[element.Nodes[i]] =
                    field.EvaluateAt(element, naturalNodes[i], standardDisplacements, enrichedDisplacements);
            }
            return nodalTensors;
        }

        private IReadOnlyDictionary<XNode2D, Tensor2D> ElementNodalTensorsExtrapolation(XContinuumElement2D element,
            double[] freeDisplacements, double[] constrainedDisplacements, IOutputField field)
        {
            throw new NotImplementedException();
        }

    }
}
