/*
 * Title:        EdgeCloudSim - Basic Edge Orchestrator implementation
 * 
 * Description: 
 * BasicEdgeOrchestrator implements basic algorithms which are
 * first/next/best/worst/random fit algorithms while assigning
 * requests to the edge devices.
 *               
 * Licence:      GPL - http://www.gnu.org/copyleft/gpl.html
 * Copyright (c) 2017, Bogazici University, Istanbul, Turkey
 */

package edu.boun.edgecloudsim.edge_orchestrator;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;


import org.apache.commons.math3.distribution.NormalDistribution;
import org.apache.commons.math3.exception.NotStrictlyPositiveException;
import org.cloudbus.cloudsim.Datacenter;
import org.cloudbus.cloudsim.core.CloudSim;

import com.sun.xml.internal.ws.policy.privateutil.PolicyUtils.Collections;

import edu.boun.edgecloudsim.core.SimManager;
import edu.boun.edgecloudsim.core.SimSettings;
import edu.boun.edgecloudsim.edge_server.EdgeHost;
import edu.boun.edgecloudsim.edge_server.EdgeVM;
import edu.boun.edgecloudsim.edge_client.CpuUtilizationModel_Custom;
import edu.boun.edgecloudsim.edge_client.Task;
import edu.boun.edgecloudsim.utils.ETCMatrix;
import edu.boun.edgecloudsim.utils.ETTMatrix;
import edu.boun.edgecloudsim.utils.Location;
import edu.boun.edgecloudsim.utils.SimLogger;

public class BasicEdgeOrchestrator extends EdgeOrchestrator {
	private int numberOfHost; //used by load balancer
	@SuppressWarnings("unused")
	private int lastSelectedHostIndex; //used by load balancer
	private int[] lastSelectedVmIndexes; //used by each host individually
	private static Datacenter receivingBS; // !!! IMPORTANT !!! DON'T USE THE METHOD Datacenter.getId(), it's messed up, use recBS instead if u need an index.
	private static int recBS = -1; //Receiving DC ID
	
	public static int getRecBS() {
		return recBS;
	}

	public void setRecBS(int recBS) {
		this.recBS = recBS;
	}

	private ArrayList<Integer> neighboringBS = new ArrayList<>(); //Neighboring BaseStations
	
	private static ETCMatrix matrix;
	private static ETTMatrix ettmatrix;
	
	public BasicEdgeOrchestrator(String _policy, String _simScenario) {
		super(_policy, _simScenario);
	}

	@Override
	public void initialize() {
		numberOfHost=SimSettings.getInstance().getNumOfEdgeHosts();
		
		lastSelectedHostIndex = -1;
		lastSelectedVmIndexes = new int[numberOfHost];
		for(int i=0; i<numberOfHost; i++)
			lastSelectedVmIndexes[i] = -1;
	}
	
	private static double getConvolveProbability(double mu, double sigma, double deadLine) {
		
		try {
		NormalDistribution distr = new NormalDistribution(mu, sigma);
		return distr.cumulativeProbability(deadLine);
		}catch(NotStrictlyPositiveException exc) {
			return 0.0;
		}
	
	}
	

	@Override
	public int getDeviceToOffload(Task task) {
		
		int result = SimSettings.GENERIC_EDGE_DEVICE_ID;
		matrix = SimLogger.getInstance().matrix;
		double dl = 0.0;
		
		//int result = SimSettings.EDGE_ORCHESTRATOR_ID;
		if(simScenario.equals("SINGLE_TIER")){
			result = SimSettings.GENERIC_EDGE_DEVICE_ID;
		}
		else {
			result = SimSettings.EDGE_ORCHESTRATOR_ID;
		}
	
		/* Calculate the location of receiving Base Station
		 * the closest one to the device(vehicle)
		 */
		
		Location deviceLoc = SimManager.getInstance().getMobilityModel().getLocation(task.getMobileDeviceId(), CloudSim.clock());
		int xdev = deviceLoc.getXPos();
		int ydev = deviceLoc.getYPos();
		
		List<Datacenter> datacenters = SimManager.getInstance().edgeServerManager.getDatacenterList();
		
		double best = 1000;

		for(int i = 0; i < datacenters.size(); i++) {
			List<EdgeHost> hostlist = datacenters.get(i).getHostList();
			
			for(EdgeHost host : hostlist) {
				Location hostLocation = host.getLocation();
				int xhost = hostLocation.getXPos();
				int yhost = hostLocation.getYPos();
				double dist = Math.sqrt((Math.pow((double)xdev-xhost, 2))+ (Math.pow((double)ydev-yhost, 2)));
				if (dist <= best) {
					best = dist;
					setReceivingBS(datacenters.get(i));
					recBS = i;
					task.setDc(i);
					dl = SimManager.getInstance().getEdgeOrchestrator().deadline(task, matrix, 0.0001,recBS);
					//task.setDeadLine(dl);
					}
				}
			}
		
		return result;
	}
	
	@Override
	public EdgeVM getVmToOffload(Task task) {
		if(simScenario.equals("TWO_TIER_WITH_EO"))
			return selectVmOnLoadBalancer(task);
		else
			return selectVmOnHost(task);
	}
	
	/*
	 * The base case policy;
	 */
	
	public EdgeVM selectVmOnHost(Task task){
		
		EdgeVM selectedVM = null;
		List<EdgeVM> vmArray = receivingBS.getVmList();
		
		for(int i = 0; i < vmArray.size(); i++) {
			double requiredCapacity = ((CpuUtilizationModel_Custom)task.getUtilizationModelCpu()).predictUtilization(vmArray.get(i).getVmType());
			double targetVmCapacity = (double)100 - vmArray.get(i).getCloudletScheduler().getTotalUtilizationOfCpu(CloudSim.clock());
			if(requiredCapacity <= targetVmCapacity)
				selectedVM = vmArray.get(i);
			}

		
		return selectedVM;
	}
	
	/*
	 * Load Balancer Policy;
	 */

	public EdgeVM selectVmOnLoadBalancer(Task task){
		EdgeVM selectedVM = null;
		
		if(policy.equalsIgnoreCase("Probability")) {
			getNeighbors(task);
			selectedVM  = getDC(task);
			
		}
		else if(policy.equalsIgnoreCase("Baseline")){
			//Location deviceLocation = SimManager.getInstance().getMobilityModel().getLocation(task.getMobileDeviceId(), CloudSim.clock());
			//int relatedHostId=deviceLocation.getServingWlanId();
			
			
			//List<EdgeVM> vmArray = SimManager.getInstance().getLocalServerManager().getVmList(relatedHostId);
			
			List<EdgeVM> vmArray =  receivingBS.getVmList();
			for(int i = 0; i < vmArray.size(); i++) {
				double requiredCapacity = ((CpuUtilizationModel_Custom)task.getUtilizationModelCpu()).predictUtilization(vmArray.get(i).getVmType());
				double targetVmCapacity = (double)100 - vmArray.get(i).getCloudletScheduler().getTotalUtilizationOfCpu(CloudSim.clock());
				if(requiredCapacity <= targetVmCapacity)
					selectedVM = vmArray.get(i);
			}
			
		}
		else if(policy.equalsIgnoreCase("MECT")) {
			
			Datacenter dcM = SimManager.getInstance().edgeServerManager.getDatacenterList().get(getMectDC(task));
			
			List<EdgeVM> vmArrayMect = dcM.getVmList();
			for(int i = 0; i < vmArrayMect.size(); i++) {
				double requiredCapacity = ((CpuUtilizationModel_Custom)task.getUtilizationModelCpu()).predictUtilization(vmArrayMect.get(i).getVmType());
				double targetVmCapacity = (double)100 - vmArrayMect.get(i).getCloudletScheduler().getTotalUtilizationOfCpu(CloudSim.clock());
				if(requiredCapacity <= targetVmCapacity) {
					selectedVM = vmArrayMect.get(i);
					}
				}
			
		}
		
		return selectedVM;
		
	}
	
	/*
	 * Calculate the deadline for the task
	 * 
	 */
	
	public double deadline(Task task, ETCMatrix b, double slack, int targetBS) {
		
		double comDelay = SimManager.getInstance().getNetworkModel().getUploadDelay(task.getMobileDeviceId(), targetBS ,task.getCloudletFileSize()) + 
							SimManager.getInstance().getNetworkModel().getDownloadDelay(targetBS, task.getMobileDeviceId(), task.getCloudletFileSize());
		double submissionTime = task.getSubmissionTime();
		double avgMu = 0;
		for(int i = 0; i < b.getDataCnum(); i++) {
			avgMu += b.getMu(i, task.getTaskType().ordinal());
		}
		double avgMuAll = avgMu/b.getDataCnum();
		double beta = 0.9;
		double alpha = 1.0;
		double deadline =  (beta*avgMuAll)+ slack + submissionTime + alpha*comDelay;
		task.setDeadLine(deadline);
		return deadline;
	}
	
	
	public double deadlineTransfer(Task task, ETCMatrix b, double slack, int targetBS, int sourceBS) {
		
		double comDelay = SimManager.getInstance().getNetworkModel().getUploadDelay(sourceBS, targetBS ,task.getCloudletFileSize()) + 
							SimManager.getInstance().getNetworkModel().getDownloadDelay(targetBS, sourceBS, task.getCloudletFileSize());
		double submissionTime = task.getSubmissionTime();
		double avgMu = 0;
		for(int i = 0; i < b.getDataCnum(); i++) {
			avgMu += b.getMu(i, task.getTaskType().ordinal());
		}
		double avgMuAll = avgMu/b.getDataCnum();
		double beta = 0.8;
		double alpha = 1.0;
		double deadline =  beta*avgMuAll+ slack + submissionTime + alpha*comDelay;
		task.setDeadLine(deadline);
		return deadline;
	}
	
	
	
	public EdgeVM getDC(Task task) {
		
		EdgeVM selectedVM = null;
		matrix = SimLogger.getInstance().matrix;
		ettmatrix = SimLogger.getInstance().ettMatrix;

		double dl1 = 0.0;
		double dl = SimManager.getInstance().getEdgeOrchestrator().deadline(task, matrix, 0.0001,recBS);
		double bestProb = matrix.getProbability(recBS, task.getTaskType().ordinal(), dl);
		
		//SimLogger.printLine(" Probability in receving DC "+ bestProb);
		
		
		
		double[] probContainer = new double[neighboringBS.size()]; 
		
		double convMu = 0.0;
		double convVariance = 0.0;
		
		for(int i = 0; i < neighboringBS.size(); i++) {
			
			convMu = matrix.getMu(neighboringBS.get(i),  task.getTaskType().ordinal())+ettmatrix.getMu(neighboringBS.get(i), task.getTaskType().ordinal());
			
			convVariance = Math.pow(matrix.getSigma(neighboringBS.get(i), task.getTaskType().ordinal()),2)+Math.pow(ettmatrix.getSigma(neighboringBS.get(i), task.getTaskType().ordinal()),2);
			
			dl1 = SimManager.getInstance().getEdgeOrchestrator().deadline(task, matrix, 0.0001,neighboringBS.get(i)); // 
			double prob; //= matrix.getProbability(neighboringBS.get(i), task.getTaskType().ordinal(), dl1);
			
			prob = getConvolveProbability(convMu, Math.sqrt(convVariance), dl1);
			
			probContainer[i] = prob;
			
				
		}
		
		double temp;
		
		for (int i = 0; i < neighboringBS.size(); i++) 
		{
			for (int j = i + 1; j < neighboringBS.size(); j++) 
			{
				if (probContainer[i] < probContainer[j]) 
				{
					temp = probContainer[i];
					probContainer[i] = probContainer[j];
					probContainer[j] = temp;
				}
			}
		}
		
		if (probContainer[0] > bestProb) {
			int dc = neighboringBS.get(0);
			List<EdgeVM> vmArray = SimManager.getInstance().edgeServerManager.getDatacenterList().get(dc).getVmList();
			
			for(int j = 0; j < vmArray.size(); j++) {
				double requiredCapacity = ((CpuUtilizationModel_Custom)task.getUtilizationModelCpu()).predictUtilization(vmArray.get(j).getVmType());
				double targetVmCapacity = (double)100 - vmArray.get(j).getCloudletScheduler().getTotalUtilizationOfCpu(CloudSim.clock());
				if(requiredCapacity <= targetVmCapacity) {
					selectedVM = vmArray.get(0);
					}
				}
		 }
		else {
			List<EdgeVM> recvmArray = receivingBS.getVmList();
			
			for(int j = 0; j < recvmArray.size(); j++) {
				double requiredCapacity = ((CpuUtilizationModel_Custom)task.getUtilizationModelCpu()).predictUtilization(recvmArray.get(j).getVmType());
				double targetVmCapacity = (double)100 - recvmArray.get(j).getCloudletScheduler().getTotalUtilizationOfCpu(CloudSim.clock());
				if(requiredCapacity <= targetVmCapacity) {
					selectedVM = recvmArray.get(j);
					}
			}
			
		}
		
		
		neighboringBS.clear();
		return selectedVM;
	}
	
	public int getMectDC(Task task) {
		
		ETCMatrix matrix = SimLogger.getInstance().matrix;
		int mectDC = recBS;
		double compValue = 99999;
		double tmpValue = 0;
		
		List<Datacenter> dcs = SimManager.getInstance().edgeServerManager.getDatacenterList();
		
		for(int i = 0; i < dcs.size(); i++) {
			tmpValue = matrix.getMu(i, task.getTaskType().ordinal());
			
			if(compValue > tmpValue) {
				compValue = tmpValue;
				mectDC = i;
				//double ddll = SimManager.getInstance().getEdgeOrchestrator().deadline(task, matrix, 0.0,mectDC);
				double mectComD = SimManager.getInstance().getNetworkModel().getUploadDelay(task.getMobileDeviceId(), mectDC ,task.getCloudletFileSize())+SimManager.getInstance().getNetworkModel().getDownloadDelay(mectDC, task.getMobileDeviceId(), task.getCloudletFileSize());
				
				//double last = ddll-mectComD;
				task.setDeadLine(task.getDeadLine()-mectComD);
				
			}
		}
		
		//SimLogger.printLine("Deadline of task "+task.getDeadLine());
		
		return mectDC;
	}
	
	/*
	 * Fill the ArrayList for neighboring Base Stations;
	 */
	
	public void getNeighbors(Task task) {
		
		List<EdgeHost> recLoc = receivingBS.getHostList();
		Location recLocation = recLoc.get(0).getLocation();
		int xRec = recLocation.getXPos();
		int yRec = recLocation.getYPos();
		
		//SimLogger.printLine(" receving bs location X :"+ xRec +" y :" +yRec);
		
		List<Datacenter> datacenters = SimManager.getInstance().edgeServerManager.getDatacenterList();

		for(int i = 0; i < datacenters.size(); i++) {
			List<EdgeHost> neighbourlist = datacenters.get(i).getHostList();
			Location neighLocation = neighbourlist.get(0).getLocation();
			int xNeigh = neighLocation.getXPos();
			int yNeigh = neighLocation.getYPos();
			
			//SimLogger.printLine(" neighbour bs location X :"+ xNeigh +" y :" +yNeigh);
			
			int xdiff = xRec-xNeigh;
			int ydiff = yRec-yNeigh;
			if(Math.abs(xdiff) == 1) {
				neighboringBS.add(i);		
			}
		}
		//SimLogger.printLine(" Neighbours of task "+task.getCloudletId() +" "+neighboringBS);
	}

	public Datacenter getReceivingBS() {
		return receivingBS;
	}

	public static void setReceivingBS(Datacenter datacenter) {
		receivingBS = datacenter;
	}
}