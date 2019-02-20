using System.Collections.Generic;
using System.Diagnostics;
using ISAAR.MSolve.FEM.Elements;
using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.LinearAlgebra.Vectors;

//TODO: should it be reused between analysis iterations (Clear method, store the node multiplicities)?
namespace ISAAR.MSolve.FEM.Postprocessing
{
    /// <summary>
    /// Recovers the nodal strains, stresses from the solution of an analysis step. For now it only works for linear analysis.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public class StrainStressField2D
    {
        private const int numTensorEntries = 3; //TODO: use a dedicated Tensor2D class

        private readonly Dictionary<Node_v2, (double[] strains, double[] stresses)> data;
        private readonly Model_v2 model;

        public StrainStressField2D(Model_v2 model)
        {
            this.model = model;
            this.data = new Dictionary<Node_v2, (double[] strains, double[] stresses)>(model.Nodes.Count);
            foreach (var node in model.Nodes)
            {
                data.Add(node, (new double[numTensorEntries], new double[numTensorEntries]));
            }
        }

        /// <summary>
        /// Calculates the strain, stress tensors at each element'sintegration points, extrapolates the its nodes and then 
        /// averages the tensors at each node over all elements it belongs to.
        /// </summary>
        /// <param name="model"></param>
        /// <param name="freeDisplacements"></param>
        /// <returns></returns>
        public void CalculateNodalTensors(IVectorView freeDisplacements)
        {
            var nodeMultiplicities = new Dictionary<Node_v2, int>();
            foreach (Node_v2 node in model.Nodes) nodeMultiplicities.Add(node, 0); // how many elements each node belongs to

            foreach (Subdomain_v2 subdomain in model.Subdomains)
            {
                foreach (Element_v2 element in subdomain.Elements)
                {
                    ContinuumElement2D elementType = (ContinuumElement2D)(element.ElementType); //TODO: remove cast

                    // Find local displacement vector
                    double[] localDisplacements = subdomain.FreeDofOrdering.ExtractVectorElementFromSubdomain(element, 
                        freeDisplacements);

                    // Calculate strains, stresses at Gauss points and extrapolate to nodes
                    (IReadOnlyList<double[]> strainsAtGPs, IReadOnlyList<double[]> stressesAtGPs) = 
                        elementType.UpdateStrainsStressesAtGaussPoints(localDisplacements);
                    IReadOnlyList<double[]> strainsAtNodes = elementType.GaussPointExtrapolation.
                        ExtrapolateTensorFromGaussPointsToNodes(strainsAtGPs, elementType.Interpolation);
                    IReadOnlyList<double[]> stressesAtNodes = elementType.GaussPointExtrapolation.
                        ExtrapolateTensorFromGaussPointsToNodes(stressesAtGPs, elementType.Interpolation);

                    // Add them to the tensors stored so far in the dictionary
                    for (int i = 0; i < elementType.Nodes.Count; ++i)
                    {
                        Node_v2 node = elementType.Nodes[i];
                        AddToTensors(node, strainsAtNodes[i], stressesAtNodes[i]);
                        ++nodeMultiplicities[node];
                    }
                }
            }

            // Divide via the node multiplicity to find the average
            foreach (Node_v2 node in model.Nodes)
            {
                int multiplicity = nodeMultiplicities[node];
                Debug.Assert(multiplicity > 0); 
                DivideTensors(node, multiplicity);
            }
        }

        public double[] GetStrainsOfNode(Node_v2 node) => data[node].strains;
        public double[] GetStressesOfNode(Node_v2 node) => data[node].stresses;
      
        private void AddToTensors(Node_v2 node, double[] strains, double[] stresses)
        {
            (double[] storedStrains, double[] storedStresses) = data[node];
            for (int i = 0; i < numTensorEntries; ++i)
            {
                storedStrains[i] += strains[i];
                storedStresses[i] += stresses[i];
            }
        }

        private void DivideTensors(Node_v2 node, double multiplicity)
        {
            (double[] storedStrains, double[] storedStresses) = data[node];
            for (int i = 0; i < numTensorEntries; ++i)
            {
                storedStrains[i] /= multiplicity;
                storedStresses[i] /= multiplicity;
            }
        }
    }
}
