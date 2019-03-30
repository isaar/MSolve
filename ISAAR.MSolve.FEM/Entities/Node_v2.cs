using ISAAR.MSolve.Discretization;
using ISAAR.MSolve.Discretization.Interfaces;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace ISAAR.MSolve.FEM.Entities
{
    public class Node_v2: INode
	{
        private readonly List<Constraint> constraints = new List<Constraint>();
        private readonly Dictionary<int, Element_v2> elementsDictionary = new Dictionary<int, Element_v2>();
        private readonly Dictionary<int, ISubdomain_v2> subdomainsDictionary = new Dictionary<int, ISubdomain_v2>();
        private readonly Dictionary<int, Subdomain_v2> nonMatchingSubdomainsDictionary = new Dictionary<int, Subdomain_v2>();

        public override string ToString()
        {
            var header = String.Format("{0}: ({1}, {2}, {3})", ID, X, Y, Z);
            string constraintsDescripton = string.Empty;
            #region removeMaria
            //UNDONE: fix text
            //foreach (var c in constraints)
            //{
            //    string con = string.Empty;
            //    switch (c)
            //    {
            //        case DOFType.Pore:
            //            con = "Pore";
            //            break;
            //        case DOFType.RotX:
            //            con = "rX";
            //            break;
            //        case DOFType.RotY:
            //            con = "rY";
            //            break;
            //        case DOFType.RotZ:
            //            con = "rZ";
            //            break;
            //        case DOFType.Unknown:
            //            con = "?";
            //            break;
            //        case DOFType.X:
            //            con = "X";
            //            break;
            //        case DOFType.Y:
            //            con = "Y";
            //            break;
            //        case DOFType.Z:
            //            con = "Z";
            //            break;
            //    }
            //    constraintsDescripton += c.ToString() + ", ";
            //}
            #endregion
            constraintsDescripton = constraintsDescripton.Length > 1 
                ? constraintsDescripton.Substring(0, constraintsDescripton.Length - 2) 
                : constraintsDescripton;

            return String.Format("{0} - Con ({1})", header, constraintsDescripton);
        }

        public List<Constraint> Constraints => constraints;

        public Dictionary<int, Element_v2> ElementsDictionary => elementsDictionary;

        public Dictionary<int, ISubdomain_v2> SubdomainsDictionary => subdomainsDictionary;

        public Dictionary<int, Subdomain_v2> NonMatchingSubdomainsDictionary => nonMatchingSubdomainsDictionary;

        //public Element EmbeddedInElement { get; set; }
        public int ID { get; set; }
        public double X { get; set; }
        public double Y { get; set; }
        public double Z { get; set; }

        public void BuildSubdomainDictionary()
        {
            foreach (Element_v2 element in elementsDictionary.Values)
                if (!subdomainsDictionary.ContainsKey(element.Subdomain.ID))
                    subdomainsDictionary.Add(element.Subdomain.ID, element.Subdomain);
        }
    }
}
