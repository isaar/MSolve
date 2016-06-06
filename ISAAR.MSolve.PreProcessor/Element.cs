using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using ISAAR.MSolve.PreProcessor.Interfaces;
using ISAAR.MSolve.Matrices.Interfaces;

namespace ISAAR.MSolve.PreProcessor
{
    public enum AbsorptionType
    {
        Unknown = 0,
        Compressional = 1,
        Shear = 2
    }

    public class Element
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

        public IList<Node> EmbeddedNodes
        {
            get { return embeddedNodes; }
        }

        public int ID { get; set; }
        public IFiniteElement ElementType { get; set; }
        //public IFiniteElementMaterial MaterialType { get; set; }
        public Subdomain Subdomain { get; set; }
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
