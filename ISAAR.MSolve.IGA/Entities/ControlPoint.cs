using ISAAR.MSolve.Discretization;
using ISAAR.MSolve.Discretization.Interfaces;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace ISAAR.MSolve.IGA.Entities
{
    public class ControlPoint : INode
	{
        protected readonly Dictionary<int, Element> elementsDictionary = new Dictionary<int, Element>();
        protected readonly Dictionary<int, Patch> patchesDictionary =new Dictionary<int, Patch>();
        private readonly List<Constraint> constraints = new List<Constraint>();

		public List<Constraint> Constrains
        {
            get { return constraints; }
        }

        public Dictionary<int, Element> ElementsDictionary
        {
            get { return elementsDictionary; }
        }

        public Dictionary<int, Patch> PatchesDictionary
        {
            get { return patchesDictionary; }
        }

        public int ID { get; set; }
        public double X { get; set; }
        public double Y { get; set; }
        public double Z { get; set; }
        public double WeightFactor { get; set; }
        public double Ksi { get; set; }
        public double Heta { get; set; }
        public double Zeta { get; set; }

        public List<Constraint> Constraints => constraints;

        public Dictionary<int, ISubdomain_v2> SubdomainsDictionary => throw new NotImplementedException();

        public void BuildPatchesDictionary()
        {
            foreach (Element element in elementsDictionary.Values)
                if (!patchesDictionary.ContainsKey(element.Patch.ID))
                    patchesDictionary.Add(element.Patch.ID, element.Patch);
        }


		public override string ToString()
        {
            var header = String.Format("{0}: ({1}, {2}, {3})", ID, X, Y, Z);
            string constrainsDescription = string.Empty;
            foreach (var c in constraints)
            {
                string con = string.Empty;
                switch (c.DOF)
                {
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
                constrainsDescription += c.ToString() + ", ";
            }
            constrainsDescription = constrainsDescription.Length>1? constrainsDescription.Substring(0, constrainsDescription.Length - 2) : constrainsDescription;

            return String.Format("{0} - Con({1})", header, constrainsDescription);
        }

		public ControlPoint Clone()
		{
			return new ControlPoint()
			{
				ID=this.ID,
				X = X,
				Y=Y,
				Z=Z,
				Ksi = Ksi,
				Heta = Heta,
				Zeta=Zeta,
				WeightFactor = WeightFactor
			};
		}

    }
}
