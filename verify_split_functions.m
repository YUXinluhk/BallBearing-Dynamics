function verify_split_functions()
% Verify split functions using data captured from OpenJones-master

% Path to captured data
dataFile = '..\OpenJones-master\verification_data.mat';
if ~exist(dataFile, 'file')
    error('Verification data file not found: %s', dataFile);
end

loaded = load(dataFile);
vdata = loaded.verification_data;

fprintf('Starting Verification...\n');
pass_count = 0;
fail_count = 0;

% --- Verify curvInner ---
if isfield(vdata, 'curvInner') && ~isempty(vdata.curvInner)
    fprintf('\nTesting curvInner (%d samples)...\n', length(vdata.curvInner));
    for i = 1:length(vdata.curvInner)
        d = vdata.curvInner(i);
        inputs = d.inputs;
        expected = d.outputs{1};

        % Call split function
        try
            actual = curvInner(inputs{1}, inputs{2}, inputs{3}, inputs{4});

            if isApproxEqual(actual, expected)
                % fprintf('Sample %d: PASS\n', i);
                pass_count = pass_count + 1;
            else
                fprintf('Sample %d: FAIL\n', i);
                fprintf('  Expected: %s\n', mat2str(expected));
                fprintf('  Actual:   %s\n', mat2str(actual));
                fail_count = fail_count + 1;
            end
        catch e
            fprintf('Sample %d: ERROR - %s\n', i, e.message);
            fail_count = fail_count + 1;
        end
    end
else
    fprintf('\nNo data for curvInner\n');
end

% --- Verify curvOuter ---
if isfield(vdata, 'curvOuter') && ~isempty(vdata.curvOuter)
    fprintf('\nTesting curvOuter (%d samples)...\n', length(vdata.curvOuter));
    for i = 1:length(vdata.curvOuter)
        d = vdata.curvOuter(i);
        inputs = d.inputs;
        expected = d.outputs{1};

        % Call split function
        try
            actual = curvOuter(inputs{1}, inputs{2}, inputs{3}, inputs{4});

            if isApproxEqual(actual, expected)
                pass_count = pass_count + 1;
            else
                fprintf('Sample %d: FAIL\n', i);
                fprintf('  Expected: %s\n', mat2str(expected));
                fprintf('  Actual:   %s\n', mat2str(actual));
                fail_count = fail_count + 1;
            end
        catch e
            fprintf('Sample %d: ERROR - %s\n', i, e.message);
            fail_count = fail_count + 1;
        end
    end
else
    fprintf('\nNo data for curvOuter\n');
end

% --- Verify calcStress ---
if isfield(vdata, 'calcStress') && ~isempty(vdata.calcStress)
    fprintf('\nTesting calcStress (%d samples)...\n', length(vdata.calcStress));
    % Mock obj for calcStress if needed, but the split function calls are static-like?
    % The split function signature: [Q, sigma_max, a, E, K, b] = calcStress(obj, delta, rho, phys)
    % obj is first arg. We can pass a dummy struct or empty.
    dummyObj = struct('cage', struct('D_p', 0)); % Minimal dummy

    for i = 1:length(vdata.calcStress)
        d = vdata.calcStress(i);
        inputs = d.inputs; % {delta, rho, phys}
        expected_outputs = d.outputs; % {Q, sigma_max, a, E, K, b}

        try
            [Q, sigma_max, a, E, K, b] = calcStress(dummyObj, inputs{1}, inputs{2}, inputs{3});

            % Check all outputs
            match = isApproxEqual(Q, expected_outputs{1}) && ...
                isApproxEqual(sigma_max, expected_outputs{2}) && ...
                isApproxEqual(a, expected_outputs{3}) && ...
                isApproxEqual(E, expected_outputs{4}) && ...
                isApproxEqual(K, expected_outputs{5}) && ...
                isApproxEqual(b, expected_outputs{6});

            if match
                pass_count = pass_count + 1;
            else
                fprintf('Sample %d: FAIL\n', i);
                if ~isApproxEqual(Q, expected_outputs{1}), fprintf('  Q Mismatch\n'); end
                if ~isApproxEqual(K, expected_outputs{5}), fprintf('  K Mismatch\n'); end
                fail_count = fail_count + 1;
            end
        catch e
            fprintf('Sample %d: ERROR - %s\n', i, e.message);
            fail_count = fail_count + 1;
        end
    end
else
    fprintf('\nNo data for calcStress\n');
end

% --- Verify cosICAngle ---
if isfield(vdata, 'cosICAngle') && ~isempty(vdata.cosICAngle)
    fprintf('\nTesting cosICAngle (%d samples)...\n', length(vdata.cosICAngle));
    for i = 1:length(vdata.cosICAngle)
        d = vdata.cosICAngle(i);
        inputs = d.inputs; % {A_2, X_2, f_i, D, delta_i}
        expected = d.outputs{1};

        try
            actual = cosICAngle(inputs{1}, inputs{2}, inputs{3}, inputs{4}, inputs{5});

            if isApproxEqual(actual, expected)
                pass_count = pass_count + 1;
            else
                fprintf('Sample %d: FAIL\n', i);
                fprintf('  Expected: %s\n', mat2str(expected));
                fprintf('  Actual:   %s\n', mat2str(actual));
                fail_count = fail_count + 1;
            end
        catch e
            fprintf('Sample %d: ERROR - %s\n', i, e.message);
            fail_count = fail_count + 1;
        end
    end
else
    fprintf('\nNo data for cosICAngle\n');
end

fprintf('\nVerification Complete.\n');
fprintf('Passed: %d\n', pass_count);
fprintf('Failed: %d\n', fail_count);
end

function res = isApproxEqual(a, b)
if isempty(a) && isempty(b)
    res = true;
    return;
end
if isempty(a) || isempty(b)
    res = false;
    return;
end

tol = 1e-10;
% Use norm if arrays
diff = abs(a - b);
res = all(diff(:) < tol);
end
