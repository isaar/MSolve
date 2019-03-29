using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.FEM.Interfaces;

//TODO: Delete this ASAP. 1) Its purpose is element-node connectivity, which should be done through interfaces and inheritance,
//      2) The order of the nodes should be defined by what is now called ElementType
namespace ISAAR.MSolve.FEM.Entities
{
    public enum AbsorptionType
    {
        Unknown = 0,
        Compressional = 1,
        Shear = 2
    }

    public class Element_v2: IElement_v2
	{
        private readonly Dictionary<int, Node_v2> nodesDictionary = new Dictionary<int, Node_v2>();
        private readonly Dictionary<DOFType, AbsorptionType> absorptions = new Dictionary<DOFType, AbsorptionType>();
        private readonly IList<Node_v2> embeddedNodes = new List<Node_v2>();

        public Dictionary<int, Node_v2> NodesDictionary => nodesDictionary;

        public Dictionary<DOFType, AbsorptionType> Absorptions => absorptions;

        IList<INode> IElement_v2.Nodes => nodesDictionary.Values.ToList<INode>();
        public IList<Node_v2> Nodes => nodesDictionary.Values.ToList<Node_v2>();


        public IList<Node_v2> EmbeddedNodes => embeddedNodes;

        public int ID { get; set; }

		IElementType_v2 IElement_v2.ElementType => ElementType;
        public IFiniteElement_v2 ElementType { get; set; }

        ISubdomain_v2 IElement_v2.Subdomain => this.Subdomain;
        public Subdomain_v2 Subdomain { get; set; }
        public int[] DOFs { get; set; }

        public void AddNode(Node_v2 node) => nodesDictionary.Add(node.ID, node);

        public void AddNodes(IList<Node_v2> nodes)
        {
            foreach (Node_v2 node in nodes) AddNode(node);
        }
    }
}
