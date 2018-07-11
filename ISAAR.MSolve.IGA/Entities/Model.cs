using ISAAR.MSolve.Numerical.LinearAlgebra.Interfaces;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.Numerical.LinearAlgebra;
using ISAAR.MSolve.IGA.Problems.Structural.Elements;
using ISAAR.MSolve.IGA.Interfaces;
using ISAAR.MSolve.IGA.Entities.Loads;
using ISAAR.MSolve.FEM.Interfaces;

namespace ISAAR.MSolve.IGA.Entities
{
    public class Model:IStructuralModel
    {
        private int totalDOFs = 0;
        private int numberOfPatches=0;
        private int numberOfInterfaces=0;
        
        private readonly Dictionary<int, ControlPoint> controlPointsDictionary = new Dictionary<int, ControlPoint>();
        private readonly Dictionary<int, Element> elementsDictionary = new Dictionary<int, Element>();
        private readonly Dictionary<int, Patch> patchesDictionary = new Dictionary<int, Patch>();

        private readonly Dictionary< int , Dictionary<DOFType,int>> controlPointsDOFsDictionary = new Dictionary<int, Dictionary<DOFType, int>>();
        private readonly IList<Load> loads = new List<Load>();
	    private readonly IList<IMassAccelerationHistoryLoad> massAccelerationHistoryLoads = new List<IMassAccelerationHistoryLoad>();
		//private readonly IList<ΙBoundaryCondition> boundaryConditions = new List<ΙBoundaryCondition>();

		#region Properties

		public Dictionary<int, ControlPoint> ControlPointsDictionary
        {
            get { return controlPointsDictionary; }
        }

        public Dictionary<int, Element> ElementsDictionary
        {
            get { return elementsDictionary; }
        }

        public Dictionary<int, Patch> PatchesDictionary
        {
            get { return patchesDictionary; }
        }

	    public Dictionary<int, ISubdomain> ISubdomainsDictionary
	    {
		    get
		    {
			    var a = new Dictionary<int, ISubdomain>();
			    foreach (var subdomain in patchesDictionary.Values)
				    a.Add(subdomain.ID, subdomain);
			    return a;
		    }
	    }

		public IList<ControlPoint> ControlPoints
        {
            get { return controlPointsDictionary.Values.ToList<ControlPoint>(); }
        }

        public IList<Element> Elements
        {
            get { return elementsDictionary.Values.ToList<Element>(); }
        }

        public IList<Patch> Patches
        {
            get { return patchesDictionary.Values.ToList<Patch>(); }
        }

        public IList<Load> Loads
        {
            get { return loads; }
        }

        //public IList<BoundaryCondition> BoundaryConditions
        //{
        //    get { return boundaryConditions; }
        //}

        public Dictionary<int, Dictionary<DOFType,int>> ControlPointDOFsDictionary
        {
            get { return controlPointsDOFsDictionary; }
        }

        public int TotalDOFs
        {
            get { return totalDOFs; }
        }

        public int NumberOfPatches
        {
            get { return numberOfPatches; }
            set { numberOfPatches = value; }
        }
        public int NumberOfInterfaces
        {
            get { return numberOfInterfaces; }
            set { numberOfInterfaces = value; }
        }

	    public IList<IMassAccelerationHistoryLoad> MassAccelerationHistoryLoads
	    {
		    get { return massAccelerationHistoryLoads; }
		}
		#endregion

		#region Data Interconnection routines

		private void BuildElementDictionaryOfEachControlPoint()
        {
            foreach (Element element in elementsDictionary.Values)
                foreach (ControlPoint controlPoint in element.ControlPoints)
                    controlPoint.ElementsDictionary.Add(element.ID, element);
        }

        private void BuildPatchOfEachElement()
        {
            foreach (Patch patch in patchesDictionary.Values)
                foreach (Element element in patch.ElementsDictionary.Values)
                    element.Patch = patch;
        }

        private void BuildInterconnectionData()
        {
            BuildPatchOfEachElement();
            BuildElementDictionaryOfEachControlPoint();
            foreach (ControlPoint controlPoint in controlPointsDictionary.Values)
                controlPoint.BuildPatchesDictionary();

            foreach (Patch patch in PatchesDictionary.Values)
            {
                patch.BuildControlPointsDictionary();
                //patch.BuildEdgesDictionary();
                //patch.BuildFacesDictionary();
            }                
        }

        private void EnumerateGlobalDOFs()
        {
            totalDOFs = 0;
            Dictionary<int, List<DOFType>> controlPointDOFTypesDictionary = new Dictionary<int, List<DOFType>>();
            foreach (Element element in elementsDictionary.Values)
            {
                for (int i = 0; i < element.ControlPoints.Count; i++)
                {
                    if (!controlPointDOFTypesDictionary.ContainsKey(element.ControlPoints[i].ID))
                        controlPointDOFTypesDictionary.Add(element.ControlPoints[i].ID, new List<DOFType>());
                    controlPointDOFTypesDictionary[element.ControlPoints[i].ID].AddRange(element.ElementType.DOFEnumerator.GetDOFTypesForDOFEnumeration(element)[i]);
                }
            }

            foreach (ControlPoint controlPoint in controlPointsDictionary.Values)
            {
                Dictionary<DOFType, int> dofsDictionary = new Dictionary<DOFType, int>();
                foreach (DOFType dofType in controlPointDOFTypesDictionary[controlPoint.ID].Distinct<DOFType>())
                {
                    int dofID = 0;
                    foreach (DOFType constraint in controlPoint.Constrains)
                    {
                        if (constraint == dofType)
                        {
                            dofID = -1;
                            break;
                        }
                    }

                    if (dofID == 0)
                    {
                        dofID = totalDOFs;
                        totalDOFs++;
                    }
                    dofsDictionary.Add(dofType, dofID);
                }
                ControlPointDOFsDictionary.Add(controlPoint.ID, dofsDictionary);
            }
        }

        private void EnumerateDOFs()
        {
            EnumerateGlobalDOFs();
            foreach (Patch patch in patchesDictionary.Values)
            {
                patch.EnumerateDOFs();
                patch.AssignGlobalControlPointDOFsFromModel(controlPointsDOFsDictionary);
            }
        }

        public void AssignLoads()
        {
            AssignControlPointLoads();
            AssignBoundaryLoads();
        }

        private void AssignBoundaryLoads()
        {
            foreach (Patch patch in patchesDictionary.Values)
            {
                foreach (Edge edge in patch.EdgesDictionary.Values)
                {
                    Dictionary<int, double> edgeLoadDictionary =edge.CalculateLoads();
                    foreach (int dof in edgeLoadDictionary.Keys)
                        if (dof!=-1) patch.Forces[dof] += edgeLoadDictionary[dof];
                }
                foreach (Face face in patch.FacesDictionary.Values)
                {
                    Dictionary<int, double> faceLoadDictionary = face.CalculateLoads();
                    foreach (int dof in faceLoadDictionary.Keys)
                        if (dof!=-1) patch.Forces[dof] += faceLoadDictionary[dof];
                }
            }
        }

        private void AssignControlPointLoads()
        {
            foreach (Patch patch in patchesDictionary.Values)
                Array.Clear(patch.Forces, 0, patch.Forces.Length);
            foreach (Load load in loads)
                foreach (Patch patch in load.ControlPoint.PatchesDictionary.Values)
                {
                    int dof = patch.ControlPointDOFsDictionary[load.ControlPoint.ID][load.DOF];
                    if (dof >= 0)
                        patch.Forces[dof] = load.Amount / load.ControlPoint.PatchesDictionary.Count;
                }
        }

        public void ConnectDataStructures()
        {
            BuildInterconnectionData();
            EnumerateDOFs();
            AssignLoads();
        }

        #endregion

        public void Clear()
        {
            loads.Clear();
            patchesDictionary.Clear();
            elementsDictionary.Clear();
            controlPointsDictionary.Clear();
            controlPointsDOFsDictionary.Clear();
        }

		public void AssignMassAccelerationHistoryLoads(int timeStep)
		{
			throw new NotImplementedException();
		}
	}
}
