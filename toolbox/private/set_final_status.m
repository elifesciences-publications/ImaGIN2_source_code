function set_final_status(status)
    if ((strcmp(status, 'OK')) || (strcmp(status, 'NOK')))
        fprintf('FINAL_STATUS: %s \n', status)
    end
end
