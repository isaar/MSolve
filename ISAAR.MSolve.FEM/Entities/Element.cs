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

        public Dictionary<int, Node> NodesDictionary
        {
            get { return nodesDictionary; }
        }

        public Dictionary<DOFType, AbsorptionType> Absorptions
        {
            get { return absorptions; }
        }

        public IList<Node> Nodes
        {
            get { return nodesDictionary.Values.ToList<Node>(); }
        }

		public IList<INode> INodes
		{
			get
			{
				IList<INode> a = new List<INode>();
				foreach (var node in nodesDictionary.Values)
					a.Add(node);
				return a;
			}
		}

		public IList<Node> EmbeddedNodes
        {
            get { return embeddedNodes; }
        }

        public int ID { get; set; }
        public IFiniteElement ElementType { get; set; }

		public IElementType IElementType
		{
			get => ElementType;
		}

        //public IFiniteElementMaterial MaterialType { get; set; }
        public Subdomain Subdomain { get; set; }
        public Subdomain_v2 Subdomain_v2 { get; set; }
        public int[] DOFs { get; set; }

        public void AddNode(Node node)
        {
            nodesDictionary.Add(node.ID, node);
        }

        public void AddNodes(IList<Node> nodes)
        {
            foreach (Node node in nodes) AddNode(node);
        }

        //public IMatrix2D<double> K
        //{
        //    get { return ElementType.StiffnessMatrix(this); }
        //}

        //public IMatrix2D<double> M
        //{
        //    get { return ElementType.MassMatrix(this); }
        //}
    }
}
