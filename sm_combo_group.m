classdef sm_combo_group < handle 
    %sm_combo_group < handle a pulsegroup that combines two different
    %groups on two sets of channels
    %   essentially, this is a wrapper around two sm_pulsegroups. 
    % it has properties:
        % groups: cell array of groups
        % n_rep: number of times each should be repeated
        % jump: where to jump to after each one
        % options: 
        % chan: channels for each pulsegroup. if not set, will get set by
        %   the consturcot
        % type: will get set by constructor. it will keep track that this
            % is an sm_combo_group and not another kind of grp
        % pack_data (transient) has wf data and markers for loading
        % last_update (transient): time the pack_data was updated
    % relevant methods:
        % constructor takes a struct with the same fields
        % compare:
        % to_struct(): returns a struct with the properties (not transient)
        % make(plsgrp,ind, clk, tbase, time)
    
    properties
        groups;
        n_rep = [];
        jump = [];
        options = '';
        chan;
    end
    
    properties (SetAccess = immutable)
      type = ''; %human readable flag, for saving/loading
    end
    
    properties (Transient=true)
        last_update;
        pack_data = struct('tab',{},'marktab',{},'wf',{},'marker',{});
    end
    
    methods
        function cg = sm_combo_group(s)
            if ~exist('s','var')|| isempty(s)
               s= struct(); 
            end
            %cg.groups = s.groups;
            if isempty(cg.type)
                cg.type = class(cg); % will populate type so we know what class made this object
            end
            fn = fieldnames(s);
            for j = 1:length(fn)
               cg.(fn{j})=s.(fn{j}); 
            end
            if isempty(cg.chan)
               for j = 1:length(cg.groups)
                  cg.chan = [cg.chan, cg.groups{j}.chan]; 
               end
            end
            cg.lint();
        end
        
        function out = compare(pg,pg2)
            if length(pg.groups)==length(pg2.groups)
               out = 1;
               for j = 1:length(pg.groups)
                  out = out && pg.groups{j}.compare(pg2.groups{j}); 
               end
            else
                out = 0;
            end
        end
        
        function b = isequal(pg,pg2)
           p1 = pg.to_struct();
           p2 = pg2.to_struct();
           b = isequal(p1,p2);
        end
        
        function pg=make(plsgrp,ind, clk, tbase, time)
           if ~exist('ind','var')
               ind = [];
           end
           if ~exist('tbase','var')
               tbase = [];
           end
           if ~exist('time','var')
               time = [];
           end
            pg = plsgrp.to_struct();
            for j = 1:length(plsgrp.groups)
               plsgrp.groups{j}.make(ind,clk,tbase,time);
            end
            plsgrp.pack_data = plsgrp.groups{1}.pack_data;
            for j = 2:length(plsgrp.groups)
                plsgrp.pack_data(j,:) = plsgrp.groups{j}.pack_data;
            end
            %plsgrp.packdata = [plsgrp.groups.pack_data];
            plsgrp.last_update = now;
        end
        
        function pg_struct = to_struct(pg)
            pg_struct = pg.groups{1}.to_struct();
            for j = 2:length(pg.groups)
                pg_struct(j) = pg.groups{j}.to_struct();
            end
        end
        
        function lint(cg)
           return; 
        end
    end %methods
end %class