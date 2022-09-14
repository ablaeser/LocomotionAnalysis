function periBout = GetBoutData( boutScans, T, loco, deform, fluor, deformVars, periParam ) % , isCSD

Nscan = numel(T);
Nbout = numel(boutScans);
NdefVars = numel(deformVars);
%runState = loco.stateBinary(:,end);
periBout = struct('Nbout',Nbout, 'dur',[], 'scan',[], 'Nscan',[], 'iso',[], 'T',[], 'Tstart',[], 'Tstop',[], 'preScan',[], 'boutScan',[], 'postScan',[], ...
    'velocity',[], 'speed',[], 'fluor',[], 'on',[], 'off',[] ); 
for v = 1:NdefVars
    periBout.(deformVars{v}) = [];
end
for b = flip(1:Nbout)
    periBout.dur(b) = T(boutScans{b}(end)) - T(boutScans{b}(1));
    % Get full bout data
    periBout.scan{b} = boutScans{b}(1)-periParam.NbaseScan:boutScans{b}(end)+periParam.NbaseScan;  % -1
    periBout.scan{b}(periBout.scan{b} < 1 | periBout.scan{b} > Nscan ) = [];
    periBout.Nscan(b) = numel(periBout.scan{b});
    periBout.T{b} = T( periBout.scan{b} ); %
    periBout.Tstart(b) = T(boutScans{b}(1));
    periBout.Tstop(b) = T(boutScans{b}(end));
    periBout.preScan{b} = find( periBout.T{b} < periBout.Tstart(b) ); % scans within bout
    periBout.boutScan{b} = find(periBout.T{b} >= periBout.Tstart(b) & periBout.T{b} <= periBout.Tstop(b)); % scans within bout
    periBout.postScan{b} = find( periBout.T{b} > periBout.Tstop(b) ); % scans within bout
    
    % Grab peri-bout locomotion, fluor, and deformation data
    periBout.velocity{b} = loco.Vdown(periBout.scan{b});
    periBout.speed{b} = loco.speedDown(periBout.scan{b});
    for v = 1:NdefVars
        periBout.(deformVars{v}){b} = deform.(deformVars{v})(periBout.scan{b},:);
    end
    periBout.fluor{b} = fluor(periBout.scan{b},:);
    
    % Concatenate onset data into arrays
    %{
    periBout.on.T = periBout.T{b} - periBout.Tstart(b); % time relative to onset of bout
    periBout.on.ind = find( periBout.on.T < periParam.on );
    periBout.on.T = periBout.on.T(periBout.on.ind); 
    periBout.on.ind = find( periBout.on.T < periParam.on );
    periBout.on.velocity(:,b) =  periBout.velocity{b}(periBout.on.ind);
    periBout.on.speed(:,b) = periBout.speed{b}(periBout.on.ind); 
    for v = 1:NdefVars
        periBout.on.(deformVars{v})(:,:,b) = periBout.(deformVars{v}){b}(periBout.on.ind,:);
    end
    periBout.on.fluor(:,:,b) = periBout.fluor{b}(periBout.on.ind,:);
    %}

    % Concatenate offset data into arrays
    %{
    periBout.off.T = periBout.T{b} - periBout.Tstop(b); % time relative to offset of running
    periBout.off.ind = find( periBout.off.T > -periParam.on ); 
    periBout.off.T = periBout.off.T(periBout.off.ind); 
    periBout.off.velocity(:,b) =  periBout.velocity{b}(periBout.off.ind);
    periBout.off.speed(:,b) = periBout.speed{b}(periBout.off.ind); 
    for v = 1:NdefVars
        periBout.off.(deformVars{v})(:,:,b) = periBout.(deformVars{v}){b}(periBout.off.ind,:);
    end
    periBout.off.fluor(:,:,b) = periBout.fluor{b}(periBout.off.ind,:);
    %}
end

for b = flip(1:Nbout)
    Trel = periBout.T{b} - periBout.Tstart(b);
    onInd = find(Trel < periParam.on );
    periBout.on.T{b} = Trel(onInd); % periBout.T{b}
end

end