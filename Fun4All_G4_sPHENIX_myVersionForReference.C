int Fun4All_G4_sPHENIX( const int nEvents = 10,
                        const float genEnergy = 5.,
                        const char* particleType = "EMinus",
                        const int iEvent = 1, 
                        const char* inputFile = "/gpfs02/phenix/prod/sPHENIX/preCDR/pro.1-beta.5/single_particle/spacal1d/fieldmap/G4Hits_sPHENIX_e-_eta0_16GeV.root",
                        const char * outputFile = "G4sPHENIXCells.root") {

    //======================
    // What to run
    //======================

    // read files in HepMC format (typically output from event generators like hijing or pythia)
    const bool readhepmc = false; // read HepMC files
    // Or:
    // Use particle generator
    const bool runpythia8 = false;

    bool do_bbc         = true;
    bool do_magnet      = true;
    bool do_preshower   = false;

    bool do_svtx        = true;
    bool do_svtx_cell   = true;
    bool do_svtx_track  = true;
    bool do_svtx_eval   = true;

    bool do_cemc        = true;
    bool do_cemc_cell   = true;
    bool do_cemc_twr    = true; // may be same as 'cell' for sPHENIX
    bool do_cemc_cluster = true;
    bool do_cemc_eval   = true;

    bool do_global          = true;
    bool do_global_fastsim  = false;

    //---------------
    // Load libraries
    //---------------
    gSystem->Load("libfun4all.so");
    gSystem->Load("libg4detectors.so");
    gSystem->Load("libphhepmc.so");
    gSystem->Load("libg4testbench.so");
    gSystem->Load("libg4hough.so");
    gSystem->Load("libcemc.so");
    gSystem->Load("libg4eval.so");
    gSystem->Load("libRawClusterBuilderIA.so");
    //gSystem->Load("libMyCaloEvaluator.so");

    // establish the geometry and reconstruction setup
    gROOT->LoadMacro("G4Setup_sPHENIX.C");
    G4Init(do_svtx, do_preshower, do_cemc, false, do_magnet, false, false);
                                        // hcalin           hcalout  pipe

    int absorberactive      = 1; // set to 1 to make all absorbers active volumes
    const string magfield   = "/phenix/upgrades/decadal/fieldmaps/sPHENIX.2d.root"; 
    const float magfield_rescale = 1.4/1.5; // scale the map to a 1.4 T field

    //---------------
    // Fun4All server
    //---------------
    Fun4AllServer *se = Fun4AllServer::instance();
    se->Verbosity(0); 

    //-----------------
    // Event generation
    //-----------------
    // Toss low multiplicity dummy events.
    PHG4SimpleEventGenerator *gen = new PHG4SimpleEventGenerator();
    if (strcmp(particleType, "EMinus") == 0)        gen->add_particles("e-", 1); 
    else if (strcmp(particleType, "Gamma") == 0)    gen->add_particles("gamma", 1); 
    else if (strcmp(particleType, "Pi0") == 0)      gen->add_particles("pi0", 1); 
    
    gen->set_vertex_distribution_function(PHG4SimpleEventGenerator::Uniform,
                                          PHG4SimpleEventGenerator::Uniform,
                                          PHG4SimpleEventGenerator::Uniform);
    gen->set_vertex_distribution_mean(0.0, 0.0, 0.0);
    gen->set_vertex_distribution_width(0.0, 0.0, 5.0);
    gen->set_vertex_size_function(PHG4SimpleEventGenerator::Uniform);
    gen->set_vertex_size_parameters(0.0, 0.0);

    // Main control of particle kinematics.
    gen->set_eta_range(0., 0.);
    gen->set_phi_range(0., 0.);
    gen->set_pt_range(genEnergy, genEnergy);
    gen->Embed(1);
    gen->Verbosity(0);
    se->registerSubsystem(gen);

    //---------------------
    // Detector description
    //---------------------
    G4Setup(absorberactive, magfield,       TPythia6Decayer::kAll,
            do_svtx,        do_preshower,   do_cemc, 
            false,          do_magnet,      false, 
            false,          magfield_rescale);

    //---------
    // BBC Reco
    //---------
    if (do_bbc) {
        gROOT->LoadMacro("G4_Bbc.C");
        BbcInit();
        Bbc_Reco();
    }

    //------------------
    // Detector Division
    //------------------
    if (do_svtx_cell) Svtx_Cells();
    if (do_cemc_cell) CEMC_Cells();

    //-----------------------------
    // CEMC towering and clustering
    //-----------------------------
    if (do_cemc_twr)     CEMC_Towers();
    //if (do_cemc_cluster) CEMC_Clusters(); // calls RawClusterBuilder().

    //--------------
    // SVTX tracking
    //--------------
    if (do_svtx_track) Svtx_Reco();

    //-----------------
    // Global Vertexing
    //-----------------
    if (do_global)  {
        gROOT->LoadMacro("G4_Global.C");
        Global_Reco();
    } else if (do_global_fastsim) {
        gROOT->LoadMacro("G4_Global.C");
        Global_FastSim();
    }  

    string particleString = particleType;
    RawClusterBuilderIA* ClusterBuilder = new RawClusterBuilderIA("EmcRawClusterBuilder");
    ClusterBuilder->Detector("CEMC");
    ClusterBuilder->SetGenPT(genEnergy);
    ClusterBuilder->SetParticleType(particleString);
    ClusterBuilder->SetEvent(iEvent);
    ClusterBuilder->ClusterSimple(false);
    se->registerSubsystem(ClusterBuilder);

    //MyCaloEvaluator* myCaloEval = new MyCaloEvaluator("CEMCEVALUATOR", "CEMC", "muhfilez.root");
    //myCaloEval->SetGenPT(genEnergy);
    //myCaloEval->SetParticleType(particleString);
    //se->registerSubsystem(myCaloEval);

    // Simulation evaluation
    //----------------------
    if (do_svtx_eval)       Svtx_Eval("g4svtx_eval.root");
    if (do_cemc_eval)       CEMC_Eval("g4cemc_eval.root");

    // Register my TowerDumper class object with Fun4AllServer. 
    //----------------------
    //TowerDumper* td = new TowerDumper();
    //se->registerSubsystem(td);

    //-------------- 
    // IO management
    //--------------
    // for single particle generators we just need something which drives
    // the event loop, the Dummy Input Mgr does just that
    Fun4AllInputManager *in = new Fun4AllDummyInputManager( "JADE");
    se->registerInputManager(in);

    //-----------------
    // Event processing
    //-----------------
    if (nEvents < 0) {
        return;
    } else if (nEvents == 0 && !readhits && !readhepmc) {
        cout << "Don't use 0 for number of events when using particle generators" << endl;
        return;
    }

    cout << "Running " << nEvents << " number of events." << endl;
    se->run(nEvents);

    //-----
    // Exit
    //-----
    se->End();
    std::cout << "All done" << std::endl;
    delete se;
    gSystem->Exit(0);
}
