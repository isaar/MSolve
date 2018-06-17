using ISAAR.MSolve.Discretization.Interfaces;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace ISAAR.MSolve.IGA.Entities
{
    public class ControlPoint:INode
	{
        private readonly List<DOFType> constrains = new List<DOFType>();
        private readonly Dictionary<int, Element> elementsDictionary = new Dictionary<int, Element>();
        private readonly Dictionary<int, Patch> patchesDictionary =new Dictionary<int, Patch>();
        
        public List<DOFType> Constrains
        {
            get { return constrains; }
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
            foreach (var c in constrains)
            {
                string con = string.Empty;
                switch (c)
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
