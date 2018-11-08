using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.FEM.Interfaces;
using ISAAR.MSolve.FEM.Elements;

namespace ISAAR.MSolve.FEM.Embedding
{
    public class Hexa8LAndNLTranslationTransformationVector : IEmbeddedDOFInHostTransformationVector
    {
        private readonly DOFType[] translationOnlyDOFTypes = new DOFType[] { DOFType.X, DOFType.Y, DOFType.Z };

        public IList<DOFType> GetDependentDOFTypes { get { return translationOnlyDOFTypes; } }

        public IList<IList<DOFType>> GetDOFTypesOfHost(EmbeddedNode node)
        {
            return node.EmbeddedInElement.ElementType.GetElementDOFTypes(node.EmbeddedInElement);
        }

        public double[][] GetTransformationVector(EmbeddedNode node)
        {
            CheckElementType(node.EmbeddedInElement.ElementType);

            const int commonDofsPerNode = 3;
            const int hostDofsPerNode = 3;
            const int hostShapeFunctionLength = 8;
            double[] hostShapeFunctions = ((IEmbeddedHostElement)node.EmbeddedInElement.ElementType).GetShapeFunctionsForNode(node.EmbeddedInElement, node);

            var transformation = new double[commonDofsPerNode][];
            for (int j = 0; j < commonDofsPerNode; j++)
            {
                transformation[j] = new double[hostShapeFunctionLength * hostDofsPerNode];
                for (int k = 0; k < hostShapeFunctionLength; k++)
                    transformation[j][hostDofsPerNode * k + j] = hostShapeFunctions[k];
            }
            
            return transformation;
        }

        private void CheckElementType(IFiniteElement element)
        {
            bool validElement = element is Hexa8;
            validElement |= element is Hexa8NonLinear;
            if (!(validElement)) throw new ArgumentException("Host element is not Hexa8 or Hexa8NL.");
        }
    }
}
