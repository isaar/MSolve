using ISAAR.MSolve.Discretization.Interfaces;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace ISAAR.MSolve.FEM.Entities
{
    public class Node: INode
	{
        private readonly List<DOFType> constraints = new List<DOFType>();
        private readonly Dictionary<int, Element> elementsDictionary = new Dictionary<int, Element>();
        private readonly Dictionary<int, Subdomain> subdomainsDictionary = new Dictionary<int, Subdomain>();
        private readonly Dictionary<int, Subdomain> nonMatchingSubdomainsDictionary = new Dictionary<int, Subdomain>();

        public override string ToString()
        {
            var header = String.Format("{0}: ({1}, {2}, {3})", ID, X, Y, Z);
            string constraintsDescripton = string.Empty;
            foreach (var c in constraints)
            {
                string con = string.Empty;
                switch (c)
                {
                    case DOFType.Pore:
                        con = "Pore";
                        break;
                    case DOFType.RotX:
                        con = "rX";
                        break;
                    case DOFType.RotY:
                        con = "rY";
                        break;
                    case DOFType.RotZ:
                        con = "rZ";
                        break;
                    case DOFType.Unknown:
                        con = "?";
                        break;
                    case DOFType.X:
                        con = "X";
                        break;
                    case DOFType.Y:
                        con = "Y";
                        break;
                    case DOFType.Z:
                        con = "Z";
                        break;
                }
                constraintsDescripton += c.ToString() + ", ";
            }
            constraintsDescripton = constraintsDescripton.Length > 1 ? constraintsDescripton.Substring(0, constraintsDescripton.Length - 2) : constraintsDescripton;

            return String.Format("{0} - Con ({1})", header, constraintsDescripton);
        }

        public List<DOFType> Constraints
        {
            get { return constraints; }
        }

        public Dictionary<int, Element> ElementsDictionary
        {
            get { return elementsDictionary; }
        }

        public Dictionary<int, Subdomain> SubdomainsDictionary
        {
            get { return subdomainsDictionary; }
        }

        public Dictionary<int, Subdomain> NonMatchingSubdomainsDictionary
        {
            get { return nonMatchingSubdomainsDictionary; }
        }

        //public Element EmbeddedInElement { get; set; }
        public int ID { get; set; }
        public double X { get; set; }
        public double Y { get; set; }
        public double Z { get; set; }

        public void BuildSubdomainDictionary()
        {
            foreach (Element element in elementsDictionary.Values)
                if (!subdomainsDictionary.ContainsKey(element.Subdomain.ID))
                    subdomainsDictionary.Add(element.Subdomain.ID, element.Subdomain);
        }
    }
}
