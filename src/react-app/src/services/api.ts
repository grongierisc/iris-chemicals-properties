const API_BASE_URL = 'http://localhost:50667';

export const uploadSDF = async (file: File) => {
    const formData = new FormData();
    formData.append('file', file);
    const response = await fetch(`${API_BASE_URL}/sdf`, {
        method: 'POST',
        body: formData,
    });
    const sessionId = response.headers.get('x-session-id');
    return { data: await response.json(), sessionId };
};

export const processSMILES = async (smiles: string) => {
    const response = await fetch(`${API_BASE_URL}/smiles`, {
        method: 'POST',
        headers: {
            'Content-Type': 'application/json',
        },
        body: JSON.stringify({ smiles }),
    });
    const sessionId = response.headers.get('x-session-id');
    return { data: await response.json(), sessionId: sessionId };
};

export const getSMILESImage = async (smiles: string) => {
    const response = await fetch(`${API_BASE_URL}/smiles/img`, {
        method: 'POST',
        headers: {
            'Content-Type': 'application/json',
        },
        body: JSON.stringify({ smiles }),
    });
    const sessionId = response.headers.get('x-session-id');
    const blob = await response.blob();
    return { data: blob, sessionId };
};

export const getSDFImage = async (file: File) => {
    const formData = new FormData();
    formData.append('file', file);
    const response = await fetch(`${API_BASE_URL}/sdf/img`, {
        method: 'POST',
        body: formData,
    });
    const sessionId = response.headers.get('x-session-id');
    const blob = await response.blob();
    return { data: blob, sessionId };
};

export const compareSMILES = async (smiles1: string, smiles2: string) => {
    const response = await fetch(`${API_BASE_URL}/compare/smiles`, {
        method: 'POST',
        headers: {
            'Content-Type': 'application/json',
        },
        body: JSON.stringify({ smiles_a: smiles1, smiles_b: smiles2 }),
    });
    const sessionId = response.headers.get('x-session-id');
    return { data: await response.json(), sessionId };
};

export const compareSMILESImage = async (smiles1: string, smiles2: string) => {
    const response = await fetch(`${API_BASE_URL}/compare/smiles/img`, {
        method: 'POST',
        headers: {
            'Content-Type': 'application/json',
        },
        body: JSON.stringify({ smiles_a: smiles1, smiles_b: smiles2 }),
    });
    const sessionId = response.headers.get('x-session-id');
    const blob = await response.blob();
    return { data: blob, sessionId };
};

export const compareSDF = async (file1: File, file2: File) => {
    const formData = new FormData();
    formData.append('file_a', file1);
    formData.append('file_b', file2);
    const response = await fetch(`${API_BASE_URL}/compare/sdf`, {
        method: 'POST',
        body: formData,
    });
    const sessionId = response.headers.get('x-session-id');
    return { data: await response.json(), sessionId };
};

export const compareSDFImage = async (file1: File, file2: File) => {
    const formData = new FormData();
    formData.append('file_a', file1);
    formData.append('file_b', file2);
    const response = await fetch(`${API_BASE_URL}/compare/sdf/img`, {
        method: 'POST',
        body: formData,
    });
    const sessionId = response.headers.get('x-session-id');
    const blob = await response.blob();
    return { data: blob, sessionId };
};

export const compareMixed = async (smiles: string, file: File) => {
    const formData = new FormData();
    formData.append('file', file);
    const response = await fetch(`${API_BASE_URL}/compare?smiles=${smiles}`, {
        method: 'POST',
        body: formData,
    });
    const sessionId = response.headers.get('x-session-id');
    return { data: await response.json(), sessionId };
};

export const compareMixedImage = async (smiles: string, file: File) => {
    const formData = new FormData();
    formData.append('file', file);
    const response = await fetch(`${API_BASE_URL}/compare/img?smiles=${smiles}`, {
        method: 'POST',
        body: formData,
    });
    const sessionId = response.headers.get('x-session-id');
    const blob = await response.blob();
    return { data: blob, sessionId };
};
