[CmdletBinding()]
param([string]$DiagnosticJson = 'results/benchmark_diagnostic_v06.json')
$ErrorActionPreference = 'Stop'
function Assert($Condition, [string]$Message) {
    if (-not $Condition) { throw $Message }
}
$snapshot = Get-Content -LiteralPath $DiagnosticJson -Raw | ConvertFrom-Json
$definitions = @(Import-Csv test_pairs.csv)
$sourceEntry = @($snapshot.manifest | Where-Object path -like 'results/*results.csv')
Assert ($sourceEntry.Count -eq 1) 'Expected one result source.'
$source = @(Import-Csv -LiteralPath $sourceEntry[0].path)
Assert ($snapshot.mode -eq 'historical_artifact_export') 'Historical origin not explicit.'
Assert ($snapshot.cases.Count -eq $definitions.Count) 'Lost input rows.'
$correct = 0
$paperCount = 0
for ($i = 0; $i -lt $definitions.Count; $i++) {
    $case = $snapshot.cases[$i]
    foreach ($field in @('host','virus','expected_status')) {
        Assert ($case.$field -ceq $definitions[$i].$field) "Changed definition row $i / $field"
    }
    Assert ($case.predicted_status -ceq $source[$i].predicted_status) 'Prediction changed.'
    Assert ($case.correct -eq ($case.predicted_status -ceq $case.expected_status)) 'Incorrect score.'
    if ($case.correct) { $correct++ }
    foreach ($field in @('host_aliases_used','literature_queries_attempted','raw_retrieved_count','candidates_before_selection')) {
        Assert ($null -eq $case.$field) "Fabricated missing historical field: $field"
    }
    Assert ($null -eq $case.cache_usage.hits) 'Invented cache hits.'
    Assert ($null -eq $case.search_status.complete) 'Invented search completeness.'
    $originalPapers = @($source[$i].paper_diagnostics | ConvertFrom-Json)
    $papers = @($case.candidate_papers_considered)
    Assert ($papers.Count -eq $originalPapers.Count) 'Lost candidate papers.'
    $paperCount += $papers.Count
    for ($j = 0; $j -lt $papers.Count; $j++) {
        # Compare all original paper diagnostics, including quotes and booleans.
        $before = $originalPapers[$j] | ConvertTo-Json -Depth 60 -Compress
        $after = $papers[$j].original_diagnostic | ConvertTo-Json -Depth 60 -Compress
        Assert ($before -ceq $after) "Altered evidence at input $i paper $j"
        Assert ($null -eq $papers[$j].evidence_scope.natural_vs_experimental) 'Invented evidence scope.'
    }
}
Assert ($correct -eq $snapshot.summary.correct) 'Summary mismatch.'
foreach ($entry in $snapshot.manifest) {
    if ($null -ne $entry.sha256) {
        # Source hashes identify the historical export, not a freeze on future
        # implementation. Continue enforcing immutable inputs/results/caches.
        if ($entry.path -like '*.py' -or $entry.path -like '*.sh' -or $entry.path -eq 'requirements.txt') {
            Assert ($entry.sha256 -match '^[0-9a-f]{64}$') "Invalid historical source digest: $($entry.path)"
        } else {
            Assert ((Get-FileHash -LiteralPath $entry.path).Hash.ToLowerInvariant() -ceq $entry.sha256) "Snapshot data changed: $($entry.path)"
        }
    }
}
$beforeHash = (Get-FileHash -LiteralPath $DiagnosticJson).Hash
$refused = $false
try { & ./export_benchmark_diagnostics.ps1 -OutputJson $DiagnosticJson }
catch {
    if ($_.Exception.Message -like 'Output already exists:*') { $refused = $true }
    else { throw }
}
Assert $refused 'Exporter did not refuse an overwrite.'
Assert ((Get-FileHash -LiteralPath $DiagnosticJson).Hash -ceq $beforeHash) 'Overwrite altered snapshot.'
Write-Output "PASS: $($definitions.Count) rows, $paperCount unaltered paper records, $correct correct; hashes, unknown fields and overwrite protection verified."
