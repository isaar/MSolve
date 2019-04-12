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

    public class Element: IElement
	{
        private readonly Dictionary<int, Node> nodesDictionary = new Dictionary<int, Node>();
        private readonly Dictionary<DOFType, AbsorptionType> absorptions = new Dictionary<DOFType, AbsorptionType>();
        private readonly IList<Node> embeddedNodes = new List<Node>();

        public Dictionary<int, Node> NodesDictionary => nodesDictionary;

        public Dictionary<DOFType, AbsorptionType> Absorptions => absorptions;

        IList<INode> IElement.Nodes => nodesDictionary.Values.ToList<INode>();
        public IList<Node> Nodes => nodesDictionary.Values.ToList<Node>();


        public IList<Node> EmbeddedNodes => embeddedNodes;

        public int ID { get; set; }

		IElementType IElement.ElementType => ElementType;
        public IFiniteElement ElementType { get; set; }

        ISubdomain IElement.Subdomain => this.Subdomain;
        public Subdomain Subdomain { get; set; }
        public int[] DOFs { get; set; }

        public void AddNode(Node node) => nodesDictionary.Add(node.ID, node);

        public void AddNodes(IList<Node> nodes)
        {
            foreach (Node node in nodes) AddNode(node);
        }
    }
}
