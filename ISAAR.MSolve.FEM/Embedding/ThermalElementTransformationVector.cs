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
    public class ThermalElementTransformationVector : IEmbeddedDOFInHostTransformationVector_v2
    {
        private readonly DOFType[] thermalDOFTypes = new DOFType[] { DOFType.Temperature };

        public IList<DOFType> GetDependentDOFTypes { get { return thermalDOFTypes; } }

        public IList<IList<DOFType>> GetDOFTypesOfHost(EmbeddedNode_v2 node)
        {
            return node.EmbeddedInElement.ElementType.GetElementDOFTypes(node.EmbeddedInElement);
        }

        public double[][] GetTransformationVector(EmbeddedNode_v2 node)
        {
            //CheckElementType(node.EmbeddedInElement.ElementType);

            const int commonDofsPerNode = 1;
            const int hostDofsPerNode = 1;
            const int hostShapeFunctionLength = 4; //TODO: Use the interpolation for this. Probably for the next line too.
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

        //TODO: delete this
        private void CheckElementType(IFiniteElement element)
        {
            bool validElement = element is Hexa8;
            validElement |= element is Hexa8NonLinear;
            if (!(validElement)) throw new ArgumentException("Host element is not Hexa8 or Hexa8NL.");
        }
    }
}

