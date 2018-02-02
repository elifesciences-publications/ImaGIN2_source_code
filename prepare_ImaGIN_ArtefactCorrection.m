function prepare_ImaGIN_ArtefactCorrection(Method, EventType, StartInterpolation, EndInterpolation, FileIn, FileOut)

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

D=ImaGIN_InterpolationFilter(S) ;

move( D, FileOut ) ;

set_final_status('OK')

end
