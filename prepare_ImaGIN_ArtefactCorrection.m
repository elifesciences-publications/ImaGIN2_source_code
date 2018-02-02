function prepare_ImaGIN_ArtefactCorrection(Method, EventType, StartInterpolation, EndInterpolation, FileIn, FileOut)

% fprintf( 1, [ 'prepare_ImaGIN_ArtefactCorrection branch artefact_correction_jd\n' ] ) ;
% fprintf( 1, [ 'Method\n' ] ) ; Method
% fprintf( 1, [ 'EventType\n' ] ) ; EventType
% fprintf( 1, [ 'StartInterpolation\n' ] ) ; StartInterpolation
% fprintf( 1, [ 'EndInterpolation\n' ] ) ; EndInterpolation
% fprintf( 1, [ 'FileIn\n' ] ) ; FileIn
% fprintf( 1, [ 'FileOut\n' ] ) ; FileOut

S.Fname=FileIn ;

S.method = Method;
S.EventType = EventType ;
if ischar(StartInterpolation)
    StartInterpolation=str2num(StartInterpolation);
end
S.StartInterpolation = StartInterpolation ;
if ischar(EndInterpolation)
    EndInterpolation=str2num(EndInterpolation);
end
S.EndInterpolation = EndInterpolation ;
% fprintf( 1, [ 'S\n' ] ) ; S
D=ImaGIN_InterpolationFilter(S) ;

move( D, FileOut ) ;

set_final_status('OK')

end
