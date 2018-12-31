using ISAAR.MSolve.IGA.Interfaces;
using ISAAR.MSolve.Numerical.LinearAlgebra.Interfaces;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.Discretization.Interfaces;

namespace ISAAR.MSolve.IGA.Entities
{
    public class Element:IElement
    {
        private readonly Dictionary<int, ControlPoint> controlPointDictionary =new Dictionary<int, ControlPoint>();

        private readonly Dictionary<int, Knot>knotsDictionary =new Dictionary<int, Knot>();

        public Dictionary<int, ControlPoint> ControlPointsDictionary
        {
            get { return controlPointDictionary; }
        }

        public IList<ControlPoint> ControlPoints
        {
            get { return controlPointDictionary.Values.ToList<ControlPoint>(); }
        }

	    public IList<INode> INodes
	    {
		    get
		    {
			    IList<INode> a = new List<INode>();
			    foreach (var controlPoint in controlPointDictionary.Values)
				    a.Add(controlPoint);
			    return a;
		    }
	    }

	    public Patch Patch { get; set; }

		public Dictionary<int, Knot> KnotsDictionary
        {
            get { return knotsDictionary; }
        }

        public IList<Knot> Knots
        {
            get { return knotsDictionary.Values.ToList<Knot>(); }
        }

        public int ID { get; set; }

        public Model Model { get; set; }

        public IIsogeometricElement ElementType { get; set; }

		public IElementType IElementType
		{
			get => ElementType;
		}

        //public Patch Patch { get; set; }
		

		public int[] DOFs { get; set; }

        public void AddControlPoint(ControlPoint controlPoint)
        {
            controlPointDictionary.Add(controlPoint.ID, controlPoint);
        }

        public void AddControlPoints(IList<ControlPoint> controlPoints)
        {
            foreach (ControlPoint controlPoint in controlPoints) AddControlPoint(controlPoint);
        }

        public void AddKnot(Knot knot)
        {
            knotsDictionary.Add(knot.ID, knot);
        }

        public void AddKnots(IList<Knot> knots)
        {
            foreach (Knot knot in knots) AddKnot(knot);
        }

    }
}
