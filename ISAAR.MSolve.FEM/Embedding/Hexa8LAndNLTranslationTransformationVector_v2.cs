using System;
using System.Collections.Generic;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.FEM.Elements;
using ISAAR.MSolve.FEM.Interfaces;

namespace ISAAR.MSolve.FEM.Embedding
{
    public class Hexa8LAndNLTranslationTransformationVector_v2 : IEmbeddedDOFInHostTransformationVector_v2
    {
        private readonly DOFType[] translationOnlyDOFTypes = new DOFType[] { DOFType.X, DOFType.Y, DOFType.Z };

        public IList<DOFType> GetDependentDOFTypes { get { return translationOnlyDOFTypes; } }

        public IList<IList<DOFType>> GetDOFTypesOfHost(EmbeddedNode_v2 node)
        {
            return node.EmbeddedInElement.ElementType.GetElementDOFTypes(node.EmbeddedInElement);
        }

        public double[][] GetTransformationVector(EmbeddedNode_v2 node)
        {
            CheckElementType(node.EmbeddedInElement.ElementType);

            const int commonDofsPerNode = 3;
            const int hostDofsPerNode = 3;
            const int hostShapeFunctionLength = 8;
            double[] hostShapeFunctions = ((IEmbeddedHostElement_v2)node.EmbeddedInElement.ElementType).GetShapeFunctionsForNode(node.EmbeddedInElement, node);

            var transformation = new double[commonDofsPerNode][];
            for (int j = 0; j < commonDofsPerNode; j++)
            {
                transformation[j] = new double[hostShapeFunctionLength * hostDofsPerNode];
                for (int k = 0; k < hostShapeFunctionLength; k++)
                    transformation[j][hostDofsPerNode * k + j] = hostShapeFunctions[k];
            }
            
            return transformation;
        }

        private void CheckElementType(IFiniteElement_v2 element)
        {
            bool validElement = element is Hexa8_v2;
            validElement |= element is Hexa8NonLinear_v2;
            if (!(validElement)) throw new ArgumentException("Host element is not Hexa8 or Hexa8NL.");
        }
    }
}
