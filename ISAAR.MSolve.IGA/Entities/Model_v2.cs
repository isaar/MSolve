using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.FEM.Interfaces;
using ISAAR.MSolve.IGA.Interfaces;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.Numerical.Commons;

//TODO: find what is going on with the dynamic loads and refactor them. That 564000000 in AssignMassAccelerationHistoryLoads()
//      cannot be correct.
namespace ISAAR.MSolve.IGA.Entities
{
	public class Model_v2 : IStructuralModel_v2
	{
		private int numberOfPatches = 0;
		private int numberOfInterfaces = 0;

		private readonly Dictionary<int, ControlPoint_v2> controlPointsDictionary = new Dictionary<int, ControlPoint_v2>();
		private readonly Dictionary<int, Element> elementsDictionary = new Dictionary<int, Element>();
		private readonly Dictionary<int, Patch_v2> patchesDictionary = new Dictionary<int, Patch_v2>();

		private readonly Dictionary<int, Dictionary<DOFType, int>> controlPointsDOFsDictionary = new Dictionary<int, Dictionary<DOFType, int>>();
		private readonly IList<Load> loads = new List<Load>();
		private readonly IList<IMassAccelerationHistoryLoad> massAccelerationHistoryLoads = new List<IMassAccelerationHistoryLoad>();

		private IGlobalFreeDofOrdering globalDofOrdering;

		//public IList<EmbeddedNode> EmbeddedNodes { get; } = new List<EmbeddedNode>();

		public IList<Cluster> Clusters => ClustersDictionary.Values.ToList();
		public Dictionary<int, Cluster> ClustersDictionary { get; } = new Dictionary<int, Cluster>();

		IReadOnlyList<IElement> IStructuralModel_v2.Elements => ElementsDictionary.Values.ToList();
		public IList<Element> Elements => ElementsDictionary.Values.ToList();
		public Dictionary<int, Element> ElementsDictionary { get; } = new Dictionary<int, Element>();
		public IList<Load> Loads { get; } = new List<Load>();
		public IList<IMassAccelerationHistoryLoad> MassAccelerationHistoryLoads { get; } = new List<IMassAccelerationHistoryLoad>();

		public IList<ControlPoint_v2> ControlPoints => controlPointsDictionary.Values.ToList();
		IReadOnlyList<INode> IStructuralModel_v2.Nodes => controlPointsDictionary.Values.ToList();
		public Dictionary<int, ControlPoint> ControlPointsDictionary { get; } = new Dictionary<int, ControlPoint>();

		IReadOnlyList<ISubdomain_v2> IStructuralModel_v2.Subdomains => patchesDictionary.Values.ToList();
		public IList<Patch_v2> Patches => patchesDictionary.Values.ToList();
		public Dictionary<int, Patch_v2> PatchesDictionary { get; } = new Dictionary<int, Patch_v2>();

		public Table<INode, DOFType, double> Constraints { get; private set; } = new Table<INode, DOFType, double>();//TODOMaria: maybe it's useless in model class

		public IGlobalFreeDofOrdering GlobalDofOrdering
		{
			get => globalDofOrdering;
			set
			{
				globalDofOrdering = value;
				foreach (Patch_v2 patch in Patches)
				{
					patch.DofOrdering = GlobalDofOrdering.SubdomainDofOrderings[patch];
					patch.Forces = Vector.CreateZero(patch.DofOrdering.NumFreeDofs);
				}
				//EnumerateSubdomainLagranges();
				//EnumerateDOFMultiplicity();
			}
		}

		public void AssignLoads()
		{
			foreach (Patch_v2 patch in PatchesDictionary.Values) patch.Forces.Clear();
			AssignControlPointLoads();
			AssignBoundaryLoads();
		}

		public void AssignMassAccelerationHistoryLoads(int timeStep)
		{
			throw new NotImplementedException();
		}

		private void AssignBoundaryLoads()
		{
			foreach (Patch_v2 patch in patchesDictionary.Values)
			{
				foreach (Edge edge in patch.EdgesDictionary.Values)
				{
					Dictionary<int, double> edgeLoadDictionary = edge.CalculateLoads();
					foreach (int dof in edgeLoadDictionary.Keys)
						if (dof != -1) patch.Forces[dof] += edgeLoadDictionary[dof];
				}
				foreach (Face face in patch.FacesDictionary.Values)
				{
					Dictionary<int, double> faceLoadDictionary = face.CalculateLoads();
					foreach (int dof in faceLoadDictionary.Keys)
						if (dof != -1) patch.Forces[dof] += faceLoadDictionary[dof];
				}
			}
		}

		private void AssignControlPointLoads()
		{
			foreach (Patch_v2 patch in PatchesDictionary.Values)
			{
				patch.ControlPointLoads = new Table<ControlPoint_v2, DOFType, double>();
			}

			foreach (Load load in Loads)
			{
				var cp = ((ControlPoint_v2) load.ControlPoint);
				double amountPerPatch = load.Amount / cp.PatchesDictionary_v2.Count;
				foreach (Patch_v2 patch in cp.PatchesDictionary_v2.Values)
				{
					bool wasNotContained = patch.ControlPointLoads.TryAdd(cp, load.DOF, amountPerPatch);
				}
			}

			//TODO: this should be done by the subdomain when the analyzer decides.
			foreach (Patch_v2 patch in PatchesDictionary.Values)
			{
				foreach ((ControlPoint node, DOFType dofType, double amount) in patch.ControlPointLoads)
				{
					int patchDofIdx = patch.DofOrdering.FreeDofs[node, dofType];
					patch.Forces[patchDofIdx] = amount;
				}
			}
		}


		//What is the purpose of this method? If someone wanted to clear the Model, they could just create a new one.
		public void Clear()
		{
			Loads.Clear();
			ClustersDictionary.Clear();
			PatchesDictionary.Clear();
			ElementsDictionary.Clear();
			ControlPointsDictionary.Clear();
			globalDofOrdering = null;
			Constraints.Clear();
			MassAccelerationHistoryLoads.Clear();
		}

		public void ConnectDataStructures()
		{
			BuildInterconnectionData();
			AssignConstraints();
		}

		//TODO: constraints should not be saved inside the nodes. As it is right now (22/11/2018) the same constraint 
		//      is saved in the node, the model constraints table and the subdomain constraints table. Furthermore,
		//      displacement control analyzer updates the subdomain constraints table only (another bad design decision).  
		//      It is too easy to access the wrong instance of the constraint. 
		private void AssignConstraints()
		{
			foreach (ControlPoint_v2 controlPoint in ControlPointsDictionary.Values)
			{
				if (controlPoint.Constraints_v2 == null) continue;
				foreach (Constraint constraint in controlPoint.Constraints_v2) Constraints[controlPoint, constraint.DOF] = constraint.Amount;
			}

			foreach (Patch_v2 patch in PatchesDictionary.Values) patch.ExtractConstraintsFromGlobal(Constraints);
		}
		
		private void BuildElementDictionaryOfEachControlPoint()
		{
			foreach (Element element in elementsDictionary.Values)
			foreach (ControlPoint controlPoint in element.ControlPoints)
				controlPoint.ElementsDictionary.Add(element.ID, element);
		}

		private void BuildInterconnectionData()//TODOMaria: maybe I have to generate the constraints dictionary for each subdomain here
		{
			BuildPatchOfEachElement();
			BuildElementDictionaryOfEachControlPoint();
			foreach (ControlPoint_v2 controlPoint in ControlPointsDictionary.Values) controlPoint.BuildPatchesDictionary_v2();
			
			foreach (Patch_v2 patch in PatchesDictionary.Values) patch.DefineControlPointsFromElements();
		}

		private void BuildPatchOfEachElement()
		{
			foreach (Patch_v2 patch in patchesDictionary.Values)
			foreach (Element element in patch.Elements)
				element.Patch_v2 = patch;
		}
		
	}
}
