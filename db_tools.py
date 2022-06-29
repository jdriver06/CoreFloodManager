
from PyQt5.QtWidgets import QDialog, QLabel, QLineEdit, QGridLayout, QPushButton
import pyodbc
from PyQt5.Qt import pyqtSignal
from PyQt5.QtCore import Qt
from enum import Enum


class DBObjects:

    class Procedures(Enum):

        SelectSurfactants = 'surfactants_schema.select_surfactants_v2'
        SelectSurfactantBlends = 'surfactants_schema.select_surfactant_blends'
        GetSurfactantStructure = 'surfactants_schema.get_surfactant_structure'
        SelectInvItems = 'surfactants_schema.select_inventory_items'
        SelectInvBlendItems = 'surfactants_schema.select_inventory_blend_items'
        SurfactantName = 'surfactants_schema.get_surfactant_name'
        SurfactantManufacturerName = 'surfactants_schema.get_surfactant_manufacturer_name'
        GetStksDils = 'surfactants_schema.get_stocks_and_dilutions'
        GetBlendStksDils = 'surfactants_schema.get_blend_stocks_and_dilutions'
        GetBlendStkInfo = 'surfactants_schema.get_blend_stock_info'
        GetBlendStkComp = 'surfactants_schema.get_blend_stock_composition'
        GetBlendInfo = 'surfactants_schema.get_blend_info'
        StockInfo = 'surfactants_schema.get_stock_info_for_pb'
        DilInfo = 'surfactants_schema.get_dilution_info_for_pb'
        LocationString = 'dbo.get_full_location_string'
        AllLocationStrings = 'dbo.get_all_full_location_strings'
        GetMTs = 'surfactants_schema.get_measurements_and_treatments'
        GetBlendMTs = 'surfactants_schema.get_blend_measurements_and_treatments'
        SetSDQty = 'surfactants_schema.set_stock_or_dilution_qty'
        SetBlendSDQty = 'surfactants_schema.set_blend_stock_or_dilution_qty'
        SetSDNotes = 'surfactants_schema.set_stock_or_dilution_notes'
        SetBlendSDNotes = 'surfactants_schema.set_blend_stock_or_dilution_notes'
        SetSDLoc = 'surfactants_schema.set_stock_or_dilution_loc'
        SetBlendSDLoc = 'surfactants_schema.set_blend_stock_or_dilution_loc'
        RemoveSurfPTag = 'surfactants_schema.remove_surfactant_project_tag'
        RemoveSurfBlendPTag = 'surfactants_schema.remove_surfactant_blend_project_tag'
        AddMTs = 'surfactants_schema.add_measurement_or_treatment'
        RemoveMT = 'surfactants_schema.remove_measurement_or_treatment'
        AddBlendMTs = 'surfactants_schema.add_blend_measurement_or_treatment'
        RemoveBlendMT = 'surfactants_schema.remove_blend_measurement_or_treatment'
        SubmitStk = 'surfactants_schema.submit_new_stock'
        UpdateQtyIfChanged = 'surfactants_schema.update_qty_if_changed'
        UpdateStkConvPerc = 'surfactants_schema.update_stock_conv_perc'
        UpdateBlendStkConvPerc = 'surfactants_schema.update_blend_stock_conv_perc'
        AddSurfactant = 'surfactants_schema.add_surfactant'
        GetUHforCompany = 'surfactants_schema.get_uncertain_hydrophobes_for_company'
        GetLowQtyStocks = 'surfactants_schema.get_low_qty_ueors_stocks'
        GetSurfProjTag = 'surfactants_schema.get_surfactant_project_tag'
        GetSurfBlendProjTag = 'surfactants_schema.get_surfactant_blend_project_tag'
        AddSurfProjTag = 'surfactants_schema.add_surfactant_project_tag'
        AddSurfBlendProjTag = 'surfactants_schema.add_surfactant_blend_project_tag'

        SearchCompanies = 'dbo.search_companies'
        SelectProjectsforCompanyName = 'dbo.select_projects_for_company_name'
        GetProjectNamefromID = 'dbo.get_project_name_from_id'
        SelectBrinesforProjectName = 'brines_schema.select_brines_for_project_name'
        SelectBrineDatabyName = 'brines_schema.select_brine_data_by_name'
        AddNewBrine = 'brines_schema.add_new_brine'
        RemoveBrine = 'brines_schema.remove_brine'
        UpdateBrine = 'brines_schema.update_brine_composition'
        AddProjectTagBrine = 'brines_schema.add_project_tag_for_brine'
        GetProjectTagsBrine = 'brines_schema.get_tagged_projects_for_brine'
        RemoveProjectTagBrine = 'brines_schema.remove_project_tag_for_brine'
        GetBrineAlias = 'brines_schema.get_brine_alias_in_project'
        UpdateBrineAlias = 'brines_schema.update_brine_alias_in_project'
        SearchBrineName = 'brines_schema.search_brine_by_name_or_alias'

        SearchChems = 'OSHA_schema.search_chemicals'
        GetGloveCodes = 'OSHA_schema.get_all_glove_codes'
        GetChemHazards = 'OSHA_schema.get_chemical_hazards_and_ppe'
        GetChemPHS = 'OSHA_schema.get_chemical_phs'
        ChemIDforName = 'OSHA_schema.get_chem_id_for_name'
        AvoidNotesSafety = 'OSHA_schema.get_avoid_notes_safety'
        ChangeAvoidNotesSafety = 'OSHA_schema.change_avoid_notes_safety'
        ChangePHS = 'OSHA_schema.change_phs'
        ChangeHazardsPPE = 'OSHA_schema.change_hazards_ppe'
        ChangeGloveCode = 'OSHA_schema.change_glove_code'
        SetSDSFPath = 'OSHA_schema.set_sds_file_path'
        GetSDSFPath = 'OSHA_schema.get_sds_file_path'
        GetPStateforName = 'OSHA_schema.get_p_state_for_chemical_name'
        CondInsertContainer = 'OSHA_schema.conditional_insert_container'
        GetContainerID = 'OSHA_schema.get_container_id'
        GetInvforChemName = 'OSHA_schema.get_inventory_for_chemical_name'
        GetEmptyInvforChemName = 'OSHA_schema.get_empty_inventory_for_chemical_name'
        AddOrderReqbyChemName = 'OSHA_schema.add_order_request_by_chemical_name'
        AddOrderReqbyChemNameFull = 'OSHA_schema.add_order_request_by_chemical_name_full'
        GetActiveReq = 'OSHA_schema.get_active_order_requests'
        GetContinersforChemName = 'OSHA_schema.get_containers_for_chemical_name'
        GetAllOrders = 'OSHA_schema.get_all_orders'
        ChemOrderRec = 'OSHA_schema.chemical_order_received'
        PlaceChemOrder = 'OSHA_schema.place_chemical_order'
        RedInvQty = 'OSHA_schema.reduce_inventory_quantity_v2'
        RedEmptyInvQty = 'OSHA_schema.reduce_empty_inventory_quantity'
        MoveInvItem = 'OSHA_schema.move_inventory_item'
        MoveEmptyInvItem = 'OSHA_schema.move_empty_inventory_item'
        GetSynsCAS = 'OSHA_schema.get_chemical_synonyms_and_cas_1_perc_1'
        GetAddCAS = 'OSHA_schema.get_additional_cas_nums'
        AddItemtoWatchlist = 'OSHA_schema.add_item_to_user_watchlist_by_name'
        GetUserWatchlist = 'OSHA_schema.get_user_watchlist'
        RemoveItemfromWatchlist = 'OSHA_schema.remove_item_from_watchlist_by_name'

        GetRooms = 'dbo.get_rooms'
        GetRoomLocs = 'dbo.get_rm_locations'
        CheckUpdateAvailable = 'dbo.is_newer_db_gui_version_available'
        CheckPMUpdateAvailable = 'dbo.is_newer_pm_gui_version_available'
        GetVerNotes = 'dbo.get_version_notes_and_required'
        GetVerNotesCurr = 'dbo.get_version_notes_and_required_after_current'
        GetPMVerNotesCurr = 'dbo.get_PM_version_notes_and_required_after_current'
        SubmitGeneralRequest = 'dbo.submit_general_request_v2'
        GetActiveGeneralReq = 'dbo.get_active_general_requests'
        GetInactiveGeneralReq = 'dbo.get_inactive_general_requests'
        SubmitReqResolution = 'dbo.submit_general_request_resolution'

        AddPolymer = 'polymers_schema.add_polymer'
        SearchPolymers = 'polymers_schema.search_polymers'
        AddPolymerInventory = 'polymers_schema.add_polymer_inventory_item'
        AddPolymerInventoryID = 'polymers_schema.add_polymer_inventory_item_by_id'
        MovePolymerInvItem = 'polymers_schema.move_polymer_inventory_item'
        GetPolymerInventory = 'polymers_schema.get_polymer_inventory'
        SetPolymerInvItemMass = 'polymers_schema.update_polymer_inventory_item_mass'

    class Params(Enum):

        SCName = 'class_name'
        IGName = 'ig_name'
        MName = 'm_name'
        SName = 'source_name'
        UHName = 'uh_name'
        HName = 'h_name'
        POMin = 'po_min'
        POMax = 'po_max'
        EOMin = 'eo_min'
        EOMax = 'eo_max'
        ACNMin = 'acn_min'
        ACNMax = 'acn_max'
        SID = 's_id'
        BlendID = 'blend_id'
        StkID = 'stock_id'
        BlendStkID = 'blend_stock_id'
        DilID = 'dil_id'
        BlendDilID = 'blend_dil_id'
        LocID = 'loc_id'
        SDID = 'sd_id'
        ISStk = 'is_stock'
        RecADate = 'rec_after_date'
        RecBDate = 'rec_before_date'
        SynADate = 'syn_after_date'
        SynBDate = 'syn_before_date'
        Qty = 'qty'
        ConvPerc = 'conv_perc'
        Notes = 'notes'
        MTTID = 'mtt_id'
        MTVal = 'val'
        MTDate = 'mt_date'
        MTID = 'mt_id'
        SurfName = 's_name'
        BlendName = 'blend_name'
        NEO = 'n_eo'
        NPO = 'n_po'
        MW = 'mw'
        NIG = 'n_ig'
        NRE = 'nre'
        POEORev = 'po_eo_rev'
        NPOEOC = 'n_poeochains'
        CompName = 'company_name'
        ProjName = 'project_name'
        ProjID = 'project_id'
        BrineName = 'brine_name'
        BrineAbbr = 'brine_abbr'
        BrineAlias = 'alias'
        Li = 'lithium'
        Na = 'sodium'
        K = 'potassium'
        Mg = 'magnesium'
        Ca = 'calcium'
        Ba = 'barium'
        Sr = 'strontium'
        Fe2 = 'ironII'
        F = 'fluoride'
        Cl = 'chloride'
        Br = 'bromide'
        HCO3 = 'bicarbonate'
        SO4 = 'sulfate'
        Na2CO3 = 'sodium_carbonate'
        NaCl = 'sodium_chloride'

        ChemName = 'chemical_name'
        UrgencyID = 'urgency_id'
        SearchStr = 'search_string'
        SDSFPath = 'sds_file_path'
        ChemID = 'chemical_id'
        UnitAbbr = 'unit_abbrev'
        Size = 'size'
        OrderRecDate = 'order_received_date'
        OrderID = 'order_id'
        ContID = 'container_id'
        OrderDate = 'order_date'
        InvID = 'inventory_id'
        EmptyInvID = 'empty_inventory_id'
        Red = 'red'
        MoveQty = 'move_qty'
        NewLocID = 'new_loc_id'
        Avoid = 'avoid'
        Safety = 'safety_details'
        IsPHS = 'is_phs'
        IsFlam = 'is_flammable'
        IsOx = 'is_oxidizer'
        IsRed = 'is_reducer'
        IsCorr = 'is_corrosive'
        IsPerox = 'is_peroxide_former'
        IsToxic = 'is_toxic'
        IsHealthHaz = 'is_health_hazard'
        IsWarning = 'is_warning'
        FaceShield = 'face_shield'
        Apron = 'apron'
        Hood = 'hood'
        GloveCodeID = 'glove_code_id'

        RoomID = 'room_id'
        v1 = 'v1'
        v2 = 'v2'
        v3 = 'v3'
        cv1 = 'cv1'
        cv2 = 'cv2'
        cv3 = 'cv3'
        RequestNotes = 'request_notes'
        GenReqID = 'general_request_id'
        ResolveNotes = 'resolve_notes'
        UserName = 'username'

        PolymerName = 'polymer_name'
        PolymerID = 'polymer_id'
        Manufacturer = 'manufacturer'
        ManufacturerName = 'manufacturer_name'
        RecCompany = 'rec_company'
        RecCompanyName = 'rec_company_name'
        PolymerType = 'polymer_type'
        AcrylateLow = 'acrylate_low'
        AcrylateHigh = 'acrylate_high'
        SulfonateLow = 'sulfonate_low'
        SulfonateHigh = 'sulfonate_high'
        MWLow = 'mw_low'
        MWHigh = 'mw_high'
        Discontinued = 'discontinued'
        IncludeDiscontinued = 'include_discontinued'
        LotNum = 'lot_num'
        RecDate = 'rec_date'
        InitMass = 'init_mass'
        CurrMass = 'current_mass'

    PPDict = {Procedures.SelectSurfactants: [Params.SurfName, Params.SCName, Params.IGName, Params.MName, Params.SName,
                                             Params.POMin, Params.POMax, Params.EOMin, Params.EOMax,
                                             Params.ACNMin, Params.ACNMax, Params.HName, Params.LocID],
              Procedures.SelectSurfactantBlends: [Params.SurfName, Params.SCName, Params.IGName, Params.MName, Params.SName,
                                             Params.POMin, Params.POMax, Params.EOMin, Params.EOMax,
                                             Params.ACNMin, Params.ACNMax, Params.HName, Params.LocID],
              Procedures.GetSurfactantStructure: [Params.SID],
              Procedures.SelectInvItems: [Params.SCName, Params.IGName, Params.MName, Params.SName,
                                          Params.POMin, Params.POMax, Params.EOMin, Params.EOMax,
                                          Params.ACNMin, Params.ACNMax, Params.HName, Params.LocID,
                                          Params.RecADate, Params.RecBDate, Params.SynADate, Params.SynBDate],
              Procedures.SelectInvBlendItems: [Params.SCName, Params.IGName, Params.MName, Params.SName,
                                               Params.POMin, Params.POMax, Params.EOMin, Params.EOMax,
                                               Params.ACNMin, Params.ACNMax, Params.HName, Params.LocID,
                                               Params.RecADate, Params.RecBDate, Params.SynADate, Params.SynBDate],
              Procedures.SurfactantName: [Params.SID],
              Procedures.SurfactantManufacturerName: [Params.SID],
              Procedures.GetStksDils: [Params.SID],
              Procedures.GetBlendStksDils: [Params.BlendID],
              Procedures.StockInfo: [Params.StkID],
              Procedures.GetBlendStkInfo: [Params.BlendStkID],
              Procedures.GetBlendStkComp: [Params.BlendStkID],
              Procedures.GetBlendInfo: [Params.BlendID],
              Procedures.DilInfo: [Params.DilID],
              Procedures.LocationString: [Params.LocID],
              Procedures.AllLocationStrings: [],
              Procedures.GetMTs: [Params.SDID, Params.ISStk],
              Procedures.GetBlendMTs: [Params.SDID, Params.ISStk],
              Procedures.SetSDQty: [Params.ISStk, Params.SDID, Params.Qty],
              Procedures.SetSDNotes: [Params.ISStk, Params.SDID, Params.Notes],
              Procedures.SetSDLoc: [Params.ISStk, Params.SDID, Params.LocID],
              Procedures.SetBlendSDQty: [Params.ISStk, Params.SDID, Params.Qty],
              Procedures.SetBlendSDNotes: [Params.ISStk, Params.SDID, Params.Notes],
              Procedures.SetBlendSDLoc: [Params.ISStk, Params.SDID, Params.LocID],
              Procedures.AddMTs: [Params.ISStk, Params.MTTID, Params.SDID, Params.MTVal, Params.MTDate],
              Procedures.RemoveMT: [Params.ISStk, Params.MTID],
              Procedures.AddBlendMTs: [Params.ISStk, Params.MTTID, Params.SDID, Params.MTVal, Params.MTDate],
              Procedures.RemoveBlendMT: [Params.ISStk, Params.MTID],
              Procedures.UpdateQtyIfChanged: [Params.ISStk, Params.SDID, Params.Qty],
              Procedures.UpdateStkConvPerc: [Params.StkID, Params.ConvPerc],
              Procedures.UpdateBlendStkConvPerc: [Params.StkID, Params.ConvPerc],
              Procedures.AddSurfactant: [Params.SurfName, Params.MName, Params.UHName, Params.HName, Params.NEO,
                                         Params.NPO, Params.IGName, Params.MW, Params.NIG, Params.NRE, Params.POEORev,
                                         Params.NPOEOC],
              Procedures.GetUHforCompany: [Params.CompName],
              Procedures.SelectProjectsforCompanyName: [Params.CompName],
              Procedures.SelectBrinesforProjectName: [Params.ProjName],
              Procedures.SelectBrineDatabyName: [Params.BrineName],
              Procedures.AddNewBrine: [Params.ProjName, Params.BrineName, Params.BrineAbbr, Params.Li, Params.Na,
                                       Params.K, Params.Mg, Params.Ca, Params.Ba, Params.Sr, Params.Fe2, Params.F,
                                       Params.Cl, Params.Br, Params.HCO3, Params.SO4, Params.Na2CO3, Params.NaCl],
              Procedures.RemoveBrine: [Params.BrineName],
              Procedures.UpdateBrine: [Params.BrineName, Params.BrineAbbr, Params.Li, Params.Na, Params.K, Params.Mg,
                                       Params.Ca, Params.Ba, Params.Sr, Params.Fe2, Params.F, Params.Cl, Params.Br,
                                       Params.HCO3, Params.SO4, Params.Na2CO3, Params.NaCl],
              Procedures.AddProjectTagBrine: [Params.BrineName, Params.ProjName, Params.BrineAlias],
              Procedures.GetProjectTagsBrine: [Params.BrineName],
              Procedures.RemoveProjectTagBrine: [Params.BrineName, Params.ProjName],
              Procedures.GetBrineAlias: [Params.BrineName, Params.ProjName],
              Procedures.UpdateBrineAlias: [Params.BrineAlias, Params.BrineName, Params.ProjName],
              Procedures.SearchBrineName: [Params.BrineName],
              Procedures.GetLowQtyStocks: [Params.Qty],
              Procedures.GetSurfProjTag: [Params.StkID],
              Procedures.GetSurfBlendProjTag: [Params.StkID],
              Procedures.RemoveSurfPTag: [Params.StkID],
              Procedures.RemoveSurfBlendPTag: [Params.BlendStkID],
              Procedures.AddSurfProjTag: [Params.StkID, Params.ProjName],
              Procedures.AddSurfBlendProjTag: [Params.BlendStkID, Params.ProjName],

              Procedures.SearchChems: [Params.SearchStr],
              Procedures.GetGloveCodes: [],
              Procedures.GetChemPHS: [Params.ChemName],
              Procedures.ChemIDforName: [Params.ChemName],
              Procedures.AvoidNotesSafety: [Params.ChemName],
              Procedures.ChangeAvoidNotesSafety: [Params.ChemName, Params.Avoid, Params.Notes, Params.Safety],
              Procedures.ChangePHS: [Params.ChemName, Params.IsPHS],
              Procedures.ChangeHazardsPPE: [Params.ChemName, Params.IsFlam, Params.IsCorr, Params.IsOx, Params.IsRed,
                                            Params.IsPerox, Params.IsToxic, Params.IsHealthHaz, Params.IsWarning,
                                            Params.FaceShield, Params.Apron, Params.Hood],
              Procedures.ChangeGloveCode: [Params.ChemName, Params.GloveCodeID],
              Procedures.SetSDSFPath: [Params.ChemName, Params.SDSFPath],
              Procedures.GetSDSFPath: [Params.ChemName],
              Procedures.GetPStateforName: [Params.ChemName],
              Procedures.CondInsertContainer: [Params.ChemID, Params.UnitAbbr, Params.Size],
              Procedures.GetContainerID: [Params.ChemID, Params.UnitAbbr, Params.Size],
              Procedures.GetInvforChemName: [Params.ChemName],
              Procedures.GetEmptyInvforChemName: [Params.ChemName],
              Procedures.AddOrderReqbyChemName: [Params.ChemName],
              Procedures.GetActiveReq: [],
              Procedures.GetContinersforChemName: [Params.ChemName],
              Procedures.GetAllOrders: [],
              Procedures.ChemOrderRec: [Params.OrderID, Params.Qty, Params.OrderRecDate, Params.LocID],
              Procedures.PlaceChemOrder: [Params.ContID, Params.Qty, Params.OrderDate],
              Procedures.RedInvQty: [Params.InvID, Params.Red],
              Procedures.RedEmptyInvQty: [Params.EmptyInvID, Params.Red],
              Procedures.MoveInvItem: [Params.InvID, Params.MoveQty, Params.NewLocID],
              Procedures.MoveEmptyInvItem: [Params.EmptyInvID, Params.MoveQty, Params.NewLocID],
              Procedures.GetSynsCAS: [Params.ChemName],
              Procedures.GetAddCAS: [Params.ChemName],
              Procedures.AddItemtoWatchlist: [Params.ChemName],
              Procedures.GetUserWatchlist: [Params.UserName],
              Procedures.RemoveItemfromWatchlist: [Params.UserName, Params.ChemName],

              Procedures.SearchCompanies: [Params.CompName],
              Procedures.GetRooms: [Params.InvID, Params.Red],
              Procedures.GetRoomLocs: [Params.RoomID],
              Procedures.CheckUpdateAvailable: [Params.v1, Params.v2, Params.v3],
              Procedures.GetVerNotes: [Params.v1, Params.v2, Params.v3],
              Procedures.SubmitGeneralRequest: [Params.RequestNotes, Params.UrgencyID],
              Procedures.GetActiveGeneralReq: [],
              Procedures.GetInactiveGeneralReq: [],
              Procedures.SubmitReqResolution: [Params.GenReqID, Params.ResolveNotes],
              Procedures.GetProjectNamefromID: [Params.ProjID],

              Procedures.AddPolymer: [Params.PolymerName, Params.Manufacturer, Params.RecCompany,
                                      Params.PolymerType, Params.AcrylateLow, Params.AcrylateHigh, Params.SulfonateLow,
                                      Params.SulfonateHigh, Params.MWLow, Params.MWHigh, Params.Discontinued,
                                      Params.Notes],
              Procedures.SearchPolymers: [Params.PolymerName, Params.ManufacturerName, Params.RecCompanyName,
                                          Params.AcrylateLow, Params.AcrylateHigh, Params.SulfonateLow,
                                          Params.SulfonateHigh, Params.MWLow, Params.MWHigh,
                                          Params.IncludeDiscontinued],
              Procedures.AddPolymerInventory: [Params.PolymerName, Params.ManufacturerName, Params.LotNum,
                                               Params.RecDate, Params.InitMass, Params.CurrMass, Params.LocID],
              Procedures.AddPolymerInventoryID: [Params.PolymerID, Params.LotNum, Params.RecDate, Params.InitMass,
                                                 Params.LocID],
              Procedures.MovePolymerInvItem: [Params.InvID, Params.LocID],
              Procedures.SetPolymerInvItemMass: [Params.InvID, Params.CurrMass],

              Procedures.GetPolymerInventory: [Params.PolymerID]}

    @staticmethod
    def are_params_in_procedure(procedure: Procedures, params: list) -> bool:

        for param in params:
            if param not in DBObjects.PPDict[procedure]:
                return False

        return True

    class Tables(Enum):

        CompaniesTable = 'dbo.companies'
        ProjectsTable = 'dbo.projects'
        SCTable = 'surfactants_schema.surfactant_class'
        IGTable = 'surfactants_schema.ionic_groups'
        UHTable = 'surfactants_schema.uncertain_hydrophobes'
        HTable = 'surfactants_schema.hydrophobes'
        RoomTable = 'dbo.rooms'
        RoomLocTable = 'dbo.rm_locations'
        SDBTable = 'surfactants_schema.surfactant_db'
        StkTable = 'surfactants_schema.surfactant_stocks'
        BlendStkTable = 'surfactants_schema.surfactant_blend_stocks'
        DilTable = 'surfactants_schema.surfactant_dilutions'
        BlendDilTable = 'surfactants_schema.surfactant_blend_dilutions'
        MTTTable = 'surfactants_schema.measurement_and_treatment_types'
        StkMTTable = 'surfactants_schema.stocks_measurements_and_treatments'
        DilMTTable = 'surfactants_schema.dilutions_measurements_and_treatments'

        ChemDB = 'OSHA_schema.chemical_db'
        SingleChemDB = 'OSHA_schema.single_component_chemical_db'
        ChemInv = 'OSHA_schema.chemical_inventory'
        InvUnits = 'OSHA_schema.inventory_units'
        ChemContainers = 'OSHA_schema.chemical_containers'
        ChemOrders = 'OSHA_schema.chemical_orders'
        ChemOrdersRec = 'OSHA_schema.chemical_orders_received'

    class Columns(Enum):

        SCID = 'sc_id'
        SCName = 'class_name'
        IGID = 'ig_id'
        IGName = 'ig_name'
        HID = 'h_id'
        HName = 'h_name'
        UHID = 'uh_id'
        UHFormula = 'h_formula'
        CTag = 'c_tag'
        SID = 's_id'
        BlendID = 'blend_id'
        StkID = 'stock_id'
        BlendStkID = 'blend_stock_id'
        DilID = 'dil_id'
        BlendDilID = 'blend_dil_id'
        Notes = 'notes'
        Qty = 'qty'
        MTTID = 'mtt_id'
        MTTName = 'mtt_name'
        MTID = 'mt_id'
        MTVal = 'val'
        MTDate = 'mt_date'
        LotN = 'lot_num'
        RecD = 'rec_date'
        SynD = 'syn_date'
        ActP = 'act_perc'
        ConvP = 'conv_perc'
        LocID = 'loc_id'
        RoomID = 'room_id'
        RoomName = 'room_name'
        LocName = 'loc_name'
        CompanyID = 'company_id'
        CompanyName = 'company_name'
        CompanyAbbr = 'company_abbrev'
        ProjectID = 'project_id'
        ProjectName = 'project_name'
        ProjectAbbr = 'company_abbrev'

        ChemName = 'chemical_name'
        CASNum = 'cas_num'
        MW = 'mw'
        PState = 'p_state'
        InvID = 'inventory_id'
        ContID = 'container_id'
        UnitAbbr = 'unit_abbrev'
        UnitName = 'unit_name'
        ChemID = 'chemical_id'
        Size = 'size'
        Ordered = 'ordered'
        OrderDate = 'order_date'
        OrderID = 'order_id'
        OrderRecDate = 'order_received_date'

    TCDict = {Tables.SCTable: [Columns.SCID, Columns.SCName],
              Tables.IGTable: [Columns.IGID, Columns.IGName, Columns.SCID],
              Tables.HTable: [Columns.HID, Columns.HName],
              Tables.UHTable: [Columns.UHID, Columns.UHFormula],
              Tables.SDBTable: [Columns.SID],
              Tables.StkTable: [Columns.SID, Columns.StkID, Columns.Notes, Columns.Qty],
              Tables.BlendStkTable: [Columns.BlendID, Columns.BlendStkID, Columns.Notes, Columns.Qty],
              Tables.DilTable: [Columns.StkID, Columns.DilID, Columns.Notes, Columns.Qty],
              Tables.BlendDilTable: [Columns.BlendStkID, Columns.BlendDilID, Columns.Notes, Columns.Qty],
              Tables.MTTTable: [Columns.MTTID, Columns.MTTName],
              Tables.StkMTTable: [Columns.MTID, Columns.MTTID, Columns.StkID, Columns.MTVal, Columns.MTDate],
              Tables.DilMTTable: [Columns.MTID, Columns.MTTID, Columns.DilID, Columns.MTVal, Columns.MTDate],
              Tables.RoomTable: [Columns.RoomID, Columns.RoomName],
              Tables.RoomLocTable: [Columns.LocID, Columns.RoomID, Columns.LocName],
              Tables.CompaniesTable: [Columns.CompanyID, Columns.CompanyName, Columns.CompanyAbbr],
              Tables.ProjectsTable: [Columns.ProjectID, Columns.CompanyID, Columns.ProjectName, Columns.ProjectAbbr],

              Tables.SingleChemDB: [Columns.CASNum, Columns.ChemName, Columns.MW, Columns.PState],
              Tables.ChemInv: [Columns.InvID, Columns.ContID, Columns.Qty, Columns.LocID],
              Tables.InvUnits: [Columns.UnitAbbr, Columns.UnitName],
              Tables.ChemContainers: [Columns.ChemID, Columns.UnitAbbr, Columns.Size, Columns.Ordered],
              Tables.ChemOrders: [Columns.ContID, Columns.Qty, Columns.OrderDate],
              Tables.ChemOrdersRec: [Columns.OrderID, Columns.Qty, Columns.OrderRecDate]}

    @staticmethod
    def is_column_in_table(column: Columns, table: Tables) -> bool:

        if column in DBObjects.TCDict[table]:
            return True

        return False


class LoginError(Exception):

    def __init__(self):
        super(LoginError, self).__init__()


def connect_to_database(server: str, database: str, username: str, password: str) -> pyodbc.Connection:

    cnxn = pyodbc.connect('DRIVER={ODBC Driver 17 for SQL Server};SERVER=' + server +
                          ';DATABASE=' + database +
                          ';UID=' + username +
                          ';PWD=' + password)

    return cnxn


def select_all_from_table(cnxn: pyodbc.Connection, table: str):

    return cnxn.cursor().execute('SELECT * from ' + table)


def select_column_from_table(cnxn: pyodbc.Connection, table: str, column: str):

    return cnxn.cursor().execute('SELECT {} from {}'.format(column, table))


def select_all_from_table_where_key_equals(cnxn: pyodbc.Connection, table: str, key: str, val: int):

    return cnxn.cursor().execute('SELECT * FROM {} WHERE {}={}'.format(table, key, val))


def select_column_from_table_where_key_equals(cnxn: pyodbc.Connection, table: str, column: str, key: str, val: int):

    return cnxn.cursor().execute('SELECT {} FROM {} WHERE {}={}'.format(column, table, key, val))


def update_table_set_column_where_key_equals(cnxn: pyodbc.Connection, table: str, column: str, set_val: str, key: str,
                                             val: int):

    qry = 'UPDATE {} SET {} = {} WHERE {} = {}'
    print(qry.format(table, column, set_val, key, val))
    return cnxn.cursor().execute('UPDATE {} SET {} = {} WHERE {} = {}'.format(table, column, set_val, key, val))


def execute_stored_procedure_with_params(cnxn: pyodbc.Connection, proc_name: str, params: dict):

    qry = 'EXECUTE ' + proc_name
    i = 0

    for key, value in params.items():

        if i > 0:
            qry += ','

        qry += ' @' + key + ' = '

        if isinstance(value, str) and value != 'NULL':
            qry += '\'' + value + '\''
        else:
            qry += str(value)

        i += 1

    # print(qry)

    return cnxn.cursor().execute(qry)


def insert_values_into_table(cnxn: pyodbc.Connection, table: str, values: list):

    qry = 'INSERT INTO ' + table + ' VALUES ('
    for value in values:
        if not isinstance(value, str):
            qry += str(value) + ', '
        else:
            qry += '\'' + value + '\', '
    qry = qry.replace('\'NULL\'', 'NULL')
    qry = qry[:-2] + ')'

    print(qry)

    cnxn.cursor().execute(qry)
    cnxn.cursor().commit()


def delete_from_table_where_key_equals(cnxn: pyodbc.Connection, table: str, key: str, value: int):

    qry = 'DELETE FROM {} WHERE {} = {}'
    print(qry.format(table, key, value))

    cnxn.cursor().execute(qry.format(table, key, value))
    cnxn.cursor().commit()


class DBConnection:

    def __init__(self, cnxn: pyodbc.Connection):
        self._cnxn = cnxn
        self._cursor = cnxn.cursor()

    def get_cursor(self):
        return self._cursor

    @staticmethod
    def check_column_and_table(column: DBObjects.Columns, table: DBObjects.Tables):

        if not DBObjects.is_column_in_table(column, table):
            raise(TypeError('The enumerated column is not in the enumerated table.'))

    def get_column_from_table(self, table: DBObjects.Tables, column: DBObjects.Columns):

        DBConnection.check_column_and_table(column, table)
        return select_column_from_table(self._cnxn, table.value, column.value).fetchall()

    def get_column_from_table_where(self, table: DBObjects.Tables, column: DBObjects.Columns,
                                    key: DBObjects.Columns, val: int):

        DBConnection.check_column_and_table(column, table)
        DBConnection.check_column_and_table(key, table)
        return select_column_from_table_where_key_equals(self._cnxn, table.value, column.value,
                                                         key.value, val).fetchall()

    def get_all_from_table_where(self, table: DBObjects.Tables, key: DBObjects.Columns, val: int):

        DBConnection.check_column_and_table(key, table)
        return select_all_from_table_where_key_equals(self._cnxn, table.value, key.value, val).fetchall()

    def set_column_where(self, table: DBObjects.Tables, column: DBObjects.Columns, set_val: str,
                         key: DBObjects.Columns, val: int):

        DBConnection.check_column_and_table(key, table)
        return update_table_set_column_where_key_equals(self._cnxn, table.value, column.value, set_val,
                                                        key.value, val).fetchall()

    def insert_values_into_table(self, table: DBObjects.Tables, *values):

        return insert_values_into_table(self._cnxn, table.value, list(values))

    def delete_from_table_where(self, table: DBObjects.Tables, key: DBObjects.Columns, value: int):

        return delete_from_table_where_key_equals(self._cnxn, table.value, key.value, value)

    def get_all_companies(self):

        return self.get_column_from_table(DBObjects.Tables.CompaniesTable, DBObjects.Columns.CompanyName)

    def exec_stored_procedure_with_params(self, procedure: DBObjects.Procedures, params: dict):

        return execute_stored_procedure_with_params(self._cnxn, procedure.value, params).fetchall()

    def exec_stored_procedure_with_params_no_fetchall(self, procedure: DBObjects.Procedures, params: dict):

        return execute_stored_procedure_with_params(self._cnxn, procedure.value, params)

    def exec_get_location_string(self, loc_id: int):

        params = {DBObjects.Params.LocID.value: loc_id}
        return self.exec_stored_procedure_with_params(DBObjects.Procedures.LocationString, params)

    def exec_get_all_location_strings(self):

        return self.exec_stored_procedure_with_params(DBObjects.Procedures.AllLocationStrings, {})

    def close(self):

        self._cnxn.close()


class LoginDialog(QDialog):

    login_signal = pyqtSignal(str, str)
    close_signal = pyqtSignal()

    def __init__(self, parent):
        super(LoginDialog, self).__init__(parent=parent)
        self.sf = parent.sf
        self.setWindowTitle('Login')
        self.setFixedSize(int(self.sf * 300), int(self.sf * 120))

        lyt = QGridLayout()

        name_txt = QLabel(parent=self, text='Name:')
        font = name_txt.font()
        font.setPointSize(12)
        name_txt.setFont(font)
        pass_txt = QLabel(parent=self, text='Password:')
        pass_txt.setFont(font)
        name_edit = QLineEdit(parent=self)
        name_edit.setFont(font)
        pass_edit = QLineEdit(parent=self)
        pass_edit.setFont(font)
        pass_edit.setEchoMode(QLineEdit.Password)
        login_button = QPushButton(parent=self, text='Login')
        login_button.setFont(font)
        login_button.clicked.connect(self.generate_login_signal)

        self.name_edit = name_edit
        self.pass_edit = pass_edit
        self.login_button = login_button

        lyt.addWidget(name_txt, 0, 0, 1, 1)
        lyt.addWidget(pass_txt, 1, 0, 1, 1)
        lyt.addWidget(name_edit, 0, 1, 1, 1)
        lyt.addWidget(pass_edit, 1, 1, 1, 1)
        lyt.addWidget(login_button, 2, 1, 1, 1)

        self.setLayout(lyt)

        self.setWindowModality(Qt.ApplicationModal)
        self.show()

    def generate_login_signal(self):

        self.login_signal.emit(self.name_edit.text(), self.pass_edit.text())

    def closeEvent(self, _):

        try:
            self.parent().login_dialog = None
        except Exception as e:
            print(e)
        finally:
            self.close_signal.emit()
